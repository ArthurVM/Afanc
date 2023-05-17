from collections import defaultdict
import json
from scipy.stats import tstd

from Afanc.utilities.runCommands import command
from Afanc.utilities.generalUtils import iupac


class Bam:
    """
    A class which takes a bam file and a variants bed, and produces a lineage profile according to the variant information within the variants bed.
    """

    def __init__(self, bam_file, prefix, threads=1):

        self.bam_file = bam_file
        self.prefix = prefix

        self.filetype = "cram" if bam_file[-5:]==".cram" else "bam"

    def bam2profile(self, ref_file, variants_bed):

        bam_profile = self._make_bam_profile(ref_file, variants_bed)

        ref_nucleotides = self._get_ref_nucleotides(ref_file, variants_bed)

        bam_profile = self._get_var_pos_cov(bam_profile, variants_bed, ref_nucleotides)

        variant_profile, variant_rows, snp_box = self._get_variant_profile(bam_profile, variants_bed)

        top_variants = self._calc_variant_profile(variant_profile, variant_rows)

        self._make_json(top_variants, ref_file, variants_bed)

        return top_variants


    def _make_bam_profile(self, ref_file, variants_bed):
        """ Generate a profile of the BAM file at positions defined by the variants bed file.

        BAM profile dictionary is of form
        {
            chr0 : { pos0 : { refN : d, sampleN : d, ... },
                     pos1 : { refN : d, sampleN : d, ... },
                     ...
                  },
            chr1 : { pos0 : { refN : d, sampleN : d, ... },
                     pos1 : { refN : d, sampleN : d, ... },
                     ...
                  },
            ...
        }

        where refN and sampleN are nucleotides from the reference and sample at the given pos,
        and d is the depth covered by that nucleotide at that position.
        """

        ## dictionary for the bam profile
        bam_profile = defaultdict(lambda : defaultdict(dict))

        ## make temporary bam file with any reads containing target regions within the variants bed file
        command(f"samtools view -Mb -L {variants_bed} {self.bam_file} -T {ref_file} > {self.prefix}.tmp.bam", "GETVARS").run_comm(0)
        command(f"samtools index {self.prefix}.tmp.bam", "GETVARS").run_comm(0)

        ## run a pileup
        runline = f"bcftools mpileup -f {ref_file} -T {variants_bed} {self.prefix}.tmp.bam -BI -a AD | bcftools call -mv | bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%GT\\t%AD]\\n'"
        pileup_stdout, pileup_stderr = command(runline, "GETVARS").run_comm_quiet(1)

        ## read and parse the reformatted bcftools output
        for l in pileup_stdout.decode().split("\n"):
            subdict = {}

            ## skip empty lines
            if l == "":
                continue

            chrom, pos, ref, alt, gt, ad = l.split("\t")
            pos = int(pos)

            alts = alt.split(",")

            ads = [int(x) for x in ad.split(",")]

            if gt == "0/0":
                subdict[ref] = ads[0]
            elif gt == "./.":
                subdict[ref] = 0
            else:
                genotypes = list([ref]+alts)

                ## capture ref allele and alt allele depths
                for i, genotype in enumerate(genotypes):
                    subdict[genotype] = ads[i]

            bam_profile[chrom][pos] = subdict

        return bam_profile


    def _get_ref_nucleotides(self, ref_file, variants_bed):
        """ Get reference nucleotides at positions defined in the variants bed
        """

        getfasta_line = f"bedtools getfasta -fi {ref_file} -bed {variants_bed}"
        getfasta_stdout, getfasta_stderr = command(getfasta_line, "GETVARS").run_comm_quiet(1)

        ## store reference variable pos nucleotides
        ref_nucleotides = {}
        for line in getfasta_stdout.decode().split("\n"):

            ## skip empty lines
            if line == "":
                continue

            if line[0] == ">":
                tmp_ = line.strip().strip(">").split(":")
                ref_chr = tmp_[0]
                ref_pos = int(tmp_[1].split("-")[1])
            else:
                ref_nucleotides[(ref_chr, ref_pos)] = line.strip().upper()

        return ref_nucleotides


    def _get_var_pos_cov(self, bam_profile, variants_bed, ref_nucleotides):
        """ Get coverage at each position defined within the variants bed file
        """

        view_line = f"samtools view -b -L {variants_bed} {self.prefix}.tmp.bam | bedtools coverage -a {variants_bed} -b - -d -sorted"
        view_stdout, view_stderr = command(view_line, "GETVARS").run_comm_quiet(1)

        variants_profile = defaultdict(list)

        for line in view_stdout.decode().split("\n"):

            ## skip empty lines
            if line == "":
                continue

            row = line.split("\t")
            chr = row[0]
            pos = int(row[2])
            var_ID = row[3]
            var_n = row[4]
            cov = int(row[-1])

            if chr in bam_profile and pos in bam_profile[chr]:
                # print(chr, pos, bam_profile[chr][pos], var_n, var_ID)

                ## skip instances where the nucleotide at this position does not match either the reference of the variant nucleotides
                if var_n not in bam_profile[chr][pos]:
                    continue

                variants_profile[var_ID].append(bam_profile[chr][pos][var_n])

            if chr not in bam_profile or pos not in bam_profile[chr]:
                bam_profile[chr][pos] = {ref_nucleotides[(chr, pos)]:cov}

        return bam_profile


    def _get_variant_profile(self, bam_profile, variants_bed):
        """ Takes the BAM profile and the variants bed and produces a dictionary of positions
        and ref/alt coverage values for variable regions defined in the variants bed.
        """

        ## open variant bed and read into a dictionary
        variant_rows = defaultdict(list)
        variant_profile = defaultdict(list)
        snp_box = []

        with open(variants_bed, 'r') as fin:
            for line in fin.readlines():

                ## skip commented lines
                if line.startswith("#"):
                    continue

                row = line.strip("\n").split("\t")
                chr, pos, var_ID, var_n = row[0], int(row[2]), row[3], row[4]

                variant_rows[var_ID].append(row)

                ## capture depths in a manner which deals with ambiguous IUPAC base designations
                ## ref_alt_tup[0] is the reference depth
                ## ref_alt_tup[1] is the summed depth of the variants which match the IUPAC base designations
                ref_alt_tup = [0, 0]

                if chr in bam_profile and pos in bam_profile[chr]:

                    ## handle ambiguous nucleotide IUPAC designations
                    for n in iupac(var_n):

                        ## check to see if the variant defining nucleotide at this variable position exists within the bam_profile dict
                        if n in bam_profile[chr][pos]:

                            ## update the ref_alt tuple with the variant depth
                            ref_alt_tup[1] += bam_profile[chr][pos][n]

                    ## infer the reference base depth and store in tuple
                    ref_alt_tup[0] = sum(list(bam_profile[chr][pos].values())) - ref_alt_tup[1]

                ## skip missing positions
                if ref_alt_tup == [0, 0]:
                    continue

                variant_profile[var_ID].append([chr, pos, row[5:], ref_alt_tup])
                snp_box.append([var_ID, chr, pos, ref_alt_tup[0], ref_alt_tup[1], ref_alt_tup[1]/sum(ref_alt_tup)])

        return variant_profile, variant_rows, snp_box


    def _calc_variant_profile(self, variant_profile, variant_rows, min_reads=5, max_std=0.15):
        """ Calculates the most likely variant profile from the variant profile
        """

        def check_variant(loc_depths_tup):
            """ checks the profile for a single variant and calculates the variant probabilities
            """
            loc_depths_box = []
            for chr, pos, description, [ref_d, alt_d] in loc_depths_tup:
                ## calculate the variant probability as the ratio between reference and alternate nucleotide depths
                ## check if the number of reads exceeds the given min threshold
                if ref_d+alt_d >= min_reads:
                    variant_prob = alt_d/(ref_d+alt_d)

                    ## add to the location depths tuple
                    loc_depths_box.append([chr, pos, description, [ref_d, alt_d], variant_prob])

            return loc_depths_box

        top_variants = defaultdict(dict)

        for var_ID, loc_depths_tup in variant_profile.items():
            ## calculate the probability of the variant allele at this position
            # variant_probs = [ alt_d/(ref_d+alt_d) for ref_d, alt_d in depths if ref_d+alt_d >= min_reads ]

            loc_depths_box = check_variant(loc_depths_tup)

            variant_probs = [v for _, _, _, _, v  in loc_depths_box]
            variant_depths = [d for _, _, _, d, _ in loc_depths_box]
            if sum(variant_probs) == 0.0:
                continue

            ## handle missing or singular datasets and calculate stdev of allelic probs
            if len(variant_probs) == 0:
                continue
            elif len(variant_probs) == 1:
                stdev = 0.0
            else:
                stdev = tstd(variant_probs)

            ## skip over instances of high standard deviation
            if stdev > max_std:
                continue

            ## get the sum of alternative and reference depths for this variant
            alt_d_sum = sum([alt_d for ref_d, alt_d in variant_depths])
            ref_d_sum = sum([ref_d for ref_d, alt_d in variant_depths])

            ## calculate the probability of the variant allele
            alt_prob = alt_d_sum/(alt_d_sum+ref_d_sum)

            if alt_prob < 0.05:
                continue

            ## capture description fields in the variants bed file
            allele_data = defaultdict(lambda : defaultdict(dict))
            for chr, pos, description, [ref_d, alt_d], variant_prob in loc_depths_box:
                allele_data[chr][pos] = {"variant_probability" : variant_prob, "reference_depth" : ref_d, "alt_depth" : alt_d, "description" : description }

            ## capture top variants
            top_variants[var_ID] = {"allele_frequency" : alt_prob, "allele_info" : allele_data}
                # print(var_ID, chr, pos, ref_d, alt_d)

        return top_variants


    def _make_json(self, top_variants, ref_file, variants_bed):
        """ Generates a JSON file containing results from variant profiling.
        """
        json_data = { "bam_file" : self.bam_file, "reference" : ref_file, "variants_bed" : variants_bed, "phylogenetics" : top_variants }

        with open(f"{self.prefix}.afanc-profile.json", "w") as fout:
            json.dump(json_data, fout, indent=4)
