import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
from hicrep.utils import readMcool
from hicrep import hicrepSCC
import numpy as np
import pandas as pd
import seaborn as sns
import argparse

import gzip
def readFile(infile):
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
        
def writeFile(outfile):
    if outfile.endswith((".gz","gzip")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--infile", required=True, help="The input file. [file1.cool,file2.cool,file3.cool,file4.cool,......]")
    parser.add_argument('-o', "--outfileprefix", required=True, help="The output file prefix.")
    parser.add_argument('-s', "--sizefile", required=True, help="chrom size file [will remove Y/chrY/chrM/chrMt/M/Mt]")
    parser.add_argument('-b', "--binsize", required=False, default=50000, type=int, help="resolution")
    parser.add_argument('-hsh', "--hsmooth", required=False, default=5, type=int, help="smooth window")
    parser.add_argument('-dm', "--dBPMax", required=False, default=2500000, type=int, help="smooth window")
    args = parser.parse_args()

    binSize = args.binsize
    h = args.hsmooth
    dBPMax = args.dBPMax
    bDownSample = False

    chromlist = []
    lengthlist = []
    fin = readFile(args.sizefile)
    for line in fin:
        tmp = line.strip().split()
        tmp[0] = tmp[0].replace("chr","")
        tmp[0] = tmp[0].replace("chr0","")
        tmp[0] = tmp[0].replace("Chr","")
        tmp[0] = tmp[0].replace("Chr0","")
        tmp[0] = tmp[0].replace("CHR","")
        tmp[0] = tmp[0].replace("CHR0","")
        chromlist.append(tmp[0])
        lengthlist.append(tmp[1])
    fin.close()
    chromlist = np.array(chromlist, dtype=str)
    lengthlist = np.array(lengthlist, dtype=int)

    files = args.infile.split(",")
    mylist = [files]
    mat = []
    for index1,in1 in enumerate(files):
        mylist.append([])
        mylist[index1+1].append(in1)
        mat.append([])
        for index2,in2 in enumerate(files):
            fmcool1 = in1
            fmcool2 = in2
            cool1, binSize1 = readMcool(fmcool1, binSize)
            cool2, binSize2 = readMcool(fmcool2, binSize)
            sccSub = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, chromlist)
            mylist[index1+1].append(sum(sccSub*lengthlist)/sum(lengthlist))
            mat[index1].append(sum(sccSub*lengthlist)/sum(lengthlist))

    fout = writeFile("{0}.cor.mat".format(args.outfileprefix))
    for ii in mylist:
        fout.write("{0}\n".format("\t".join([str(i) for i in ii])))
    fout.close()

    matrix = pd.DataFrame(mat,columns=files)
    matrix.index = files
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap="coolwarm", vmin=0)
    plt.savefig("{0}.pdf".format(args.outfileprefix), dpi=150, format="pdf")