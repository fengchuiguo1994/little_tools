import sys

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

import pysam
def readSam(insamfile):
    if insamfile.endswith(".bam"):
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

infile = sys.argv[1]
outfile = sys.argv[2]
fin = readSam(infile)
fout = writeFile(outfile)

myhash = {}
for read in fin:
    if read.is_unmapped:
        if -10 not in myhash:
            myhash[-10] = {}
        if -10 not in myhash[-10]:
            myhash[-10][-10] = 0
        myhash[-10][-10] += 1
        continue
    if read.mapping_quality not in myhash:
        myhash[read.mapping_quality] = {}
    diff = read.get_tag("AS") - read.get_tag("XS")
    if diff not in myhash[read.mapping_quality]:
        myhash[read.mapping_quality][diff] = 0
    myhash[read.mapping_quality][diff] += 1

for q in sorted(myhash.keys()):
    for d in sorted(myhash[q]):
        fout.write("{0}\t{1}\t{2}\n".format(q,d,myhash[q][d]))
fout.close()
fin.close()