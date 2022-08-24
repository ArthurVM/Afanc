import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import path


def modifyFasta(fasta_file, accession):
    """ Takes a fasta file and modifies each chromosome header to include
    the file name for easy grepping out of a SAM file.
    """
    records = []

    if fasta_file.endswith(".gz"):
        fin = gzip.open(fasta_file, "rt")
    else:
        fin = open(fasta_file, "r")

    for rec in SeqIO.parse(fin, "fasta"):
        record = SeqRecord(
            Seq(rec.seq),
            id=f"{rec.id}___{accession}",
            name=rec.id,
            description=rec.description,
        )
        records.append(record)

    SeqIO.write(records, f"{accession}.tmp.fa", "fasta")
    fin.close()
