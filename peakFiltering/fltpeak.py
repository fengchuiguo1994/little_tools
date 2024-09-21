import pyBigWig
import sys
import numpy as np

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

bw = pyBigWig.open(sys.argv[1])
chrom_dict = bw.chroms()

fin = readFile(sys.argv[2])
fout1 = writeFile("{0}.drop.txt".format(sys.argv[2]))
fout2 = writeFile("{0}.start.txt".format(sys.argv[2]))
fout3 = writeFile("{0}.end.txt".format(sys.argv[2]))
fout4 = writeFile("{0}.ok.txt".format(sys.argv[2]))
for line in fin:
    tmp = line.strip().split()
    if tmp[0] == "chrM":
        continue
    aa = bw.values(tmp[0],int(tmp[1]),int(tmp[2]))

    start = None
    end = None
    for i in range(len(aa)-1):
        if aa[i] == 0:
            continue
        if aa[i+1] >= aa[i]*5:
            start = i
            break
    if start != None:
        for i in range(len(aa)-1,start,-1):
            if aa[i] == 0:
                continue
            if aa[i-1] > aa[i]*5:
                end = i
                break
    else:
        for i in range(len(aa)-1,0,-1):
            if aa[i] == 0:
                continue
            if aa[i-1] > aa[i]*5:
                end = i
                break
    if start != None and end != None:
        fout1.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(line.strip(), start, aa[start], aa[start+1], end, aa[end-1], aa[end]))
    elif start == None and end == None:
        fout4.write(line)
    elif start != None and end == None:
        fout2.write("{0}\t{1}\t{2}\t{3}\n".format(line.strip(), start, aa[start], aa[start+1]))
    elif start == None and end != None:
        fout3.write("{0}\t{1}\t{2}\t{3}\n".format(line.strip(), end, aa[end-1], aa[end]))
    else:
        sys.exit(line)
fin.close()
fout1.close()
fout2.close()
fout3.close()
fout4.close()