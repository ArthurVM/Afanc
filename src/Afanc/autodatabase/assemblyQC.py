import pandas as pd
import numpy as np
import math
import re
from shutil import move
from collections import defaultdict
from scipy.stats import mode
from os import path, mkdir, chdir, listdir

from Afanc.utilities.runCommands import command


def mash(args, taxon_id, fastas):
    """ run mash
    """

    mashdist_out = path.abspath(f"{taxon_id}_mashdist.txt")
    mash_sketchline = f"mash sketch -o ref {' '.join(fastas)}"
    mash_distline = f"mash dist ref.msh ref.msh > {mashdist_out}"
    command(mash_sketchline, "MASH").run_comm_quiet(0, args.stdout, args.stderr)
    command(mash_distline, "MASH").run_comm_quiet(0, args.stdout, args.stderr)

    return mashdist_out


def buildMatrix(args, inMash):

    # read mash results and create array
    df = pd.read_csv(inMash, header=None, sep='\t')
    iniArray = df.to_numpy()

    # get tax ID
    regex = re.compile(r'\d+')
    tax = regex.findall(str(inMash))

    # calculate number of fastas, and use to split iniArray and create mash matrix
    numFastas = math.sqrt(len(iniArray))
    colList = np.split(iniArray[:,2], numFastas)
    mashMatrix = np.vstack(colList)

    # find average distance of each row
    sumArray = mashMatrix.sum(axis=1)
    avArray = sumArray/numFastas

    # get filenames and append to average distances
    filenames = iniArray[0:int(numFastas), 0]
    fileMash = np.column_stack((filenames, mashMatrix, avArray))
    fileavArray = np.column_stack((filenames, avArray))

    # save mash matrix with average distance appended on the end
    mashList = str(tax[0]) + "_mash.txt"
    np.savetxt(mashList, fileMash, delimiter=" ", fmt="%s")

    # if there are fewer than three samples in taxon then skip quality control
    if len(avArray) < 3:
        warning_txt = str(tax[0]) + "_warning.txt"
        with open(warning_txt, "w") as fout:
            print(f"Warning: Taxon {tax[0]} has fewer than 3 samples. Samples have been added to database without quality control.", file=fout)

        cleanList = str(tax[0]) + "_clean.txt"

        for elem in filenames:
            move(elem, args.cleanFasta_WDir)

        return None, None, None

    # remove duplicate files
    calcArray = fileavArray[np.unique(fileavArray[:, 1], return_index=True)[1]]

    # round avArray to 2 s.f
    rdArray = []
    for elem in avArray:
        elemrd = round(elem, -int(np.floor(np.sign(elem) * np.log10(abs(elem)))) + 2)
        rdArray.append(elemrd)

    # find the mode
    modeArray = mode(rdArray)
    modeVal = modeArray[0]

    return calcArray, tax, modeVal


def fastaMove(args, calcArray, tax, modeVal, modeRange):

    # get range around the mode
    percentRange = float(modeRange) * modeVal
    low = modeVal - percentRange
    high = modeVal + percentRange

    # collect clean fastas
    cleanIndices = np.argwhere((calcArray[:,1] >= low) & (calcArray[:,1] <= high))

    cleanFasta = []
    # if cleanIndices is not empty
    if cleanIndices.size:
        for elem in np.nditer(cleanIndices):
            cleanFasta.append(calcArray[elem, 0])

    # write to file list of high quality assemblies
    cleanList = str(tax[0]) + "_clean.txt"
    with open(cleanList, "w") as file_out:
        for elem in cleanFasta:
            move(elem, args.cleanFasta_WDir)
            # file_out.write("%s\n" % elem)
