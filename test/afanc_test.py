""" Testing module for Afanc.

Takes a database constructed using the Afanc autodatabase module,
and a csv file containing paired fastq paths and a truth value of form:

/path/to/fastq/test1_1.fq,/path/to/fastq/test1_2.fq,Mycobaterium tuberculosis
/path/to/fastq/test2_1.fq,/path/to/fastq/test2_2.fq,Mycobaterium avium
/path/to/fastq/test3_1.fq,/path/to/fastq/test3_2.fq,Mycobaterium tuberculosis Lineage 4.1

The results of the screen are compared to the truth value, and given a
designation for their accuracy:

PASS : Only 1 hit is reported and it matches the truth value.
TYPE-1 INNACURACY : There are multiple hits but the top one matches the truth value.
TYPE-2 INNACURACY : There are multiple hits and one of the non-top hits matches the truth value.
TYPE-3 INNACURACY : There is only 1 hit and the species matches the species of the truth value, indicating variant level ID failure.
TYPE-4 INNACURACY : There are multiple hits and the species of the top hit matches the species of truth value, indicating a variant level ID failure.
FAIL : Failure to successfully identify the sample based on the truth value.
RUNFAIL : The Afanc run terminated with a non-0 exit code.
"""

import os, sys, subprocess
import json
import argparse
from argparse import RawTextHelpFormatter

def run_test(args):

    ## tsv for results
    tsv_handle = open("./Afanc_test.tsv", 'w')
    print(f"Sample ID\tTruth Value\tResults\tClassification", file=tsv_handle)

    ## collect results
    results_dict = { "PASS" : 0,
                     "TYPE-1 INNACURACY" : 0,
                     "TYPE-2 INNACURACY" : 0,
                     "TYPE-3 INNACURACY" : 0,
                     "TYPE-4 INNACURACY" : 0,
                     "FAIL" : 0,
                     "RUNFAIL" : 0 }

    with open(args.run_list, 'r') as fin:
        samples_box = [ line.strip("\n").split(",") for line in fin.readlines() ]

    for fq1, fq2, truth_val in samples_box:

        if len(truth_val.split(" ")) > 2:
            variant_flag = True
        else:
            variant_flag = False

        if fq1.endswith("fq.gz"):
            prefix = os.path.basename(fq1).split("_1.fq.gz")[0]
        elif fq1.endswith("fastq.gz"):
            prefix = os.path.basename(fq1).split("_1.fastq.gz")[0]
        else:
            print("FASTQ files must end with either .fq.gz or .fastq.gz. Exiting.")
            sys.exit(1)

        runline = f"afanc screen -o {prefix} {args.db} {fq1} {fq2} -c > {prefix}.log"
        # print(runline)

        if not os.path.exists(prefix):
            subprocess.call(runline, shell=True)

        results_json = f"./{prefix}/{prefix}.json"

        if os.path.exists(results_json):
            with open(results_json, "r") as fin:
                jdata = json.load(fin)

                result, hits = check_hits(jdata["results"]["Detection_events"], truth_val)

        else:
            result = "RUNFAIL"
            print(f"{prefix} :: TRUTH={truth_val}\t{result}\n")
            ## increment the results dictionary
            results_dict[result]+=1
            ## update the results tsv
            print(f"{prefix}\t{truth_val}\tNA\t{result}", file=tsv_handle)
            continue

        ## increment the results dictionary
        results_dict[result]+=1

        print(f"{prefix} :: TRUTH={truth_val}\tRESULT={hits}\t{result}\n")
        print(f"{prefix}\t{truth_val}\t{hits}\t{result}", file=tsv_handle)

    for classification, count in results_dict.items():
        print(f"{classification} : {count}")

    tsv_handle.close()

def check_hits(hits, truth_val):
    """ checks hits against a known truth value
    """

    hit_box = []

    for hit in hits:
        if "closest_variant" in hit:
            hit = hit["closest_variant"]
        hit_box.append([hit["name"], hit["percentage"]])

    hit_box = sorted(hit_box, key=lambda x: x[1], reverse=True)

    names_box = [name for (name, perc) in hit_box]

    if len(hit_box) == 1:
        ## return a PASS if there is only 1 hit and it is the truth value
        if hit_box[0][0] == truth_val:
            return "PASS", hit_box[0][0]
        ## return a T3 INNACURACY if there is only 1 hit and the species matches the species of the truth value
        elif " ".join(hit_box[0][0].split(" ")[:2]) == " ".join(truth_val.split(" ")[:2]):
            return "TYPE-3 INNACURACY", hit_box[0][0]
        ## return a FAIL if there is only 1 hit and it does not match the truth value or its broader species
        else:
            return "FAIL", hit_box[0][0]

    elif len(hit_box) > 1:
        ## return a T1 INNACURACY if there are multiple hits but the top one matches the truth value
        if hit_box[0][0] == truth_val:
            return "TYPE-1 INNACURACY", ",".join(names_box)
        ## return a T4 INNACURACY if there are multiple hits but the species of the top one matches the species of truth value
        elif " ".join(hit_box[0][0].split(" ")[:2]) == " ".join(truth_val.split(" ")[:2]):
            return "TYPE-4 INNACURACY", ",".join(names_box)
        ## return a T2 INNACURACY if there are multiple hits and one of the non-top hits matches the truth value
        elif truth_val in names_box:
            return "TYPE-2 INNACURACY", ",".join(names_box)
        ## return a FAIL if there are multiple hits and none match the truth value
        ## N.B. this does not check each hit for species level matches
        else:
            return "FAIL", ",".join(names_box)


def parse_args(argv):
    """ argument parser for Afanc tester
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument('script', type=os.path.abspath, action='store', help=argparse.SUPPRESS)
    parser.add_argument("db", action='store', help="Database for testing.")
    parser.add_argument("run_list", action='store', help="A tab separated txt file containing the paired end .fq paths and a truth value to compare the output to.")

    run_test(parser.parse_args(argv))

parse_args(sys.argv)
