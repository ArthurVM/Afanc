import os
import sys
from collections import defaultdict


def readAccFile(db_accessions):
    dbdict = defaultdict(list)
    with open(db_accessions, "r") as fin:
        for line in fin.readlines():
            taxid, acc, name = line.strip("\n").split("\t")
            dbdict[taxid] = [acc, name]

    return dbdict


def getAccessions(report, dbdict):

    mycohits = []

    with open(report, 'r') as fin:
        for line in fin.readlines():
            if "taxon_id" in line:
                ncbitaxID = line.split('\"')[-2]
                if ncbitaxID in dbdict:
                    mycohits.append(dbdict[ncbitaxID][0])

    return mycohits
