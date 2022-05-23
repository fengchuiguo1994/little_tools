import sys
from math import log

import gzip
def readFile(infile):
    """
    infile: input file
    return: file handle
    """
    if infile.endswith((".gz","gzip",".GZ",".GZIP")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin

def writeFile(outfile):
    """
    outfile: output file
    return: file handle
    """
    if outfile.endswith((".gz","gzip",".GZ",".GZIP")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

def calcShannonEnt(dataSet):
    shannonEnt = 0.0
    total = sum(dataSet)
    for k in dataSet:
        prob = k/total
        if prob != 0.0:
            shannonEnt -= prob * log(prob,2)
    return shannonEnt

fin = readFile(sys.argv[1])
fout = writeFile(sys.argv[2])

for line in fin:
    tmp = line.strip().split()
    dataSet = []
    aa = float(tmp[1])
    flagSame = True
    for i in tmp[1:]:
        dataSet.append(float(i))
        if aa != dataSet[-1]:
            flagSame = False
    if flagSame:
        for i in range(len(dataSet)):
            dataSet[i] = 1.0
        fout.write("{0}\t{1:.4f}\n".format(line.strip(),log(len(dataSet),2) - calcShannonEnt(dataSet)))
        # fout.write("{0}\t{1}\n".format(line.strip(),log(len(dataSet),2) - calcShannonEnt(dataSet)))
    else:
        fout.write("{0}\t{1:.4f}\n".format(line.strip(),log(len(dataSet),2) - calcShannonEnt(dataSet)))
        # fout.write("{0}\t{1}\n".format(line.strip(),log(len(dataSet),2) - calcShannonEnt(dataSet)))
    
fin.close()
fout.close()