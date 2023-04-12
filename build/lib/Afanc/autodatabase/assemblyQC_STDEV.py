import math
import re
from shutil import move
from collections import defaultdict
from os import path, mkdir, chdir, listdir

import pandas as pd
import numpy as np
from scipy.stats import mode, tstd
from scipy import mean

from Afanc.utilities.runCommands import command


def mash(args, taxon_id, fastas):
    """ run mash
    """

    mashdist_out = path.abspath(f"{taxon_id}_mashdist.txt")

    ## take only the first 250 fastas for mash distance
    if len(fastas) > 250:
        print(f"More than 250 genome sequences provided at taxon {taxon_id}. Only taking the first 250 for database construction. Please manually curate the sequences within this taxon.")
        fastas = fastas[:250]

    ## run mash
    mash_sketchline = f"mash sketch -o ref {' '.join(fastas)}"
    mash_distline = f"mash dist ref.msh ref.msh > {mashdist_out}"
    command(mash_sketchline, "MASH").run_comm_quiet(0, args.stdout, args.stderr)
    command(mash_distline, "MASH").run_comm_quiet(0, args.stdout, args.stderr)

    return mashdist_out


def buildMatrix(args, inMash, average_type="mode"):

    ## read mash results and create array
    df = pd.read_csv(inMash, header=None, sep='\t')
    iniArray = df.to_numpy()

    ## get tax ID
    regex = re.compile(r'\d+')
    tax = regex.findall(str(inMash))

    ## calculate number of fastas, and use to split iniArray and create mash matrix
    numFastas = math.sqrt(len(iniArray))
    colList = np.split(iniArray[:,2], numFastas)
    mashMatrix = np.vstack(colList)

    ## find average distance of each row
    sumArray = mashMatrix.sum(axis=1)
    avArray = sumArray/numFastas

    ## get filenames and append to average distances
    filenames = iniArray[0:int(numFastas), 0]
    fileMash = np.column_stack((filenames, mashMatrix, avArray))
    fileavArray = np.column_stack((filenames, avArray))

    ## save mash matrix with average distance appended on the end
    mashList = str(tax[0]) + "_mash.txt"
    np.savetxt(mashList, fileMash, delimiter=" ", fmt="%s")

    ## if there are fewer than three samples in taxon then skip quality control
    if len(avArray) < 3:
        warning_txt = str(tax[0]) + "_warning.txt"
        with open(warning_txt, "w") as fout:
            print(f"Warning: Taxon {tax[0]} has fewer than 3 samples. Samples have been added to database without quality control.", file=fout)

        cleanList = str(tax[0]) + "_clean.txt"

        return None, None, None

    ## remove duplicate files
    # calcArray = fileavArray[np.unique(fileavArray[:, 1], return_index=True)[1]]
    calcArray = fileavArray

    ## round avArray to 2 s.f
    rdArray = []
    for elem in avArray:
        elemrd = round(elem, -int(np.floor(np.sign(elem) * np.log10(abs(elem)))) + 2)
        rdArray.append(elemrd)

    ## find the centroid (mean or mode)
    if average_type == "mean":
        centroid = np.array([mean(rdArray)])
    elif average_type == "mode":
        centroid = mode(rdArray)[0]

    return calcArray, tax, centroid


def fastaMove(args, calcArray, tax, centroid, capture_range):
    """ Move high quality assemblies to the cleanFasta directory
    """

    ## get user defined number of SDs around the centroid
    # perc_range = tstd(calcArray[:,1])*capture_range
    perc_range = float(capture_range) * centroid
    low = centroid - perc_range
    high = centroid + perc_range

    # collect clean fastas
    cleanIndices = np.argwhere((calcArray[:,1] >= low) & (calcArray[:,1] <= high))

    fasta_box = []
    # if cleanIndices is not empty
    if cleanIndices.size:
        for elem in np.nditer(cleanIndices):
            fasta_box.append(calcArray[elem, 0])

    # write to file list of high quality assemblies
    clean_fastas = []
    cleanList = str(tax[0]) + "_clean.txt"
    with open(cleanList, "w") as file_out:
        for fasta in fasta_box:
            move(fasta, args.cleanFasta_WDir)
            clean_fastas.append(path.join(args.cleanFasta_WDir, fasta.split("/")[-1]))
            file_out.write(f"{fasta}\n")

    return clean_fastas
