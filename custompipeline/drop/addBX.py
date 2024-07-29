import sys

### Operating common files  ###
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

class FqRead2:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip()
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()

fin1 = FqRead2(sys.argv[1])
fin2 = FqRead2(sys.argv[2])
outprefix = sys.argv[3]
fout1 = writeFile("{0}.BX.read1.fastq".format(outprefix))
fout2 = writeFile("{0}.BX.read2.fastq".format(outprefix))

for ((rid1,rseq1,rsyb1,rqual1),(rid2,rseq2,rsyb2,rqual2)) in zip(fin1,fin2):
    tmp1 = rid1.split("_")
    tmp2 = rid2.split("_")
    readID = tmp1[0]
    bxid = tmp1[1]
    readid = "{0}-{1}".format(readID,bxid)
    rseq1 = rseq1[33:] # 0:33 Tn5
    rqual1 = rqual1[33:]
    fout1.write("{0}\n{1}\n{2}\n{3}\n".format(readid,rseq1,rsyb1,rqual1))
    fout2.write("{0}\n{1}\n{2}\n{3}\n".format(readid,rseq2,rsyb2,rqual2))