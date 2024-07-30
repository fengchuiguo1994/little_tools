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

fin = readSam(sys.argv[1])
# fout = writeFile(sys.argv[2])
fout = writeSam(sys.argv[2], header=fin.header)
fdrop = writeSam("{0}.drop.bam".format(sys.argv[2]), header=fin.header)
if len(sys.argv) > 3:
    dis = int(sys.argv[3])
else:
    dis = 8000

n = 0
m = 0
total = 0
for read1 in fin:
    read2 = next(fin)
    total += 1
    if read1.reference_name == read2.reference_name:
        if abs(read1.isize) != abs(read2.isize):
            print(read1)
            print(read2)
            sys.exit("error")
        if abs(read1.isize) < 8000:
            fdrop.write(read1)
            fdrop.write(read2)
            n += 1
        else:
            fout.write(read1)
            fout.write(read2)
            m += 1
    else:
        fout.write(read1)
        fout.write(read2)
        m += 1
fout.close()
fdrop.close()
print(n)
print(m)
print(total)