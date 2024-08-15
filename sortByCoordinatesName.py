import sys
from collections import OrderedDict

import pysam
def readSam(insamfile):
    if insamfile.endswith(".bam"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam.gz"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam"):
        insam = pysam.AlignmentFile(insamfile,'r')
    else:
        raise ValueError("the input sam/bam file is not end with sam or bam!")
    return insam
        
def writeSam(outsamfile,header):
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

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

fin = readSam(sys.argv[1])
mydict = {}
for read in fin:
    if read.query_name not in mydict:
        mydict[read.query_name] = 0
    mydict[read.query_name] += 1
fin.close()

fin = readSam(sys.argv[1])
fout = writeSam(sys.argv[2], header=fin.header)
readdict = {}
for read in fin:
    if read.query_name not in readdict:
        readdict[read.query_name] = [read]
    else:
        readdict[read.query_name].append(read)
        if len(readdict[read.query_name]) == mydict[read.query_name]:
            for k in readdict[read.query_name]:
                fout.write(k)
            del(readdict[read.query_name])
            del(mydict[read.query_name])
fin.close()
fout.close()
# fout = writeFile(sys.argv[2])
# for rid in mydict.keys():
#     fout.write("{0}\t{1}\n".format(rid, mydict[rid]))
# fout.close()