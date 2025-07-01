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
class FqRead:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip().split()[0]
            rid = rid[1:]
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()

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

myset = set()
fin = readSam(sys.argv[1])
for read in fin.fetch("chr19"):
    if read.reference_name == "chr19":
        if read.mapping_quality > 10:
            myset.add(read.query_name)
            if len(myset) == 1000000:
                break
fin.close()

fin = FqRead(sys.argv[2])
fout = writeFile(sys.argv[3])
for rid,rseq,rsyb,rqual in fin:
    if rid in myset:
        fout.write("@{0}\n{1}\n{2}\n{3}\n".format(rid,rseq,rsyb,rqual))
fout.close()