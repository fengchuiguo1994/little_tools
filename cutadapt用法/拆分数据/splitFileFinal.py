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
        fout22 = gzip.open(outfile,'wt')
    else:
        fout22 = open(outfile,'w')
    return fout22

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

class linkerinfoRead:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        flag = None
        outlist = []
        for line in self.fin:
            if flag == None:
                flag = line.strip().split()[0]
            elif flag != line.strip().split()[0]:
                yield flag,outlist[0]
                flag = line.strip().split()[0]
                outlist = []
            outlist.append(line.strip())
        self.fin.close()
        yield flag,outlist[0]

info1 = linkerinfoRead(sys.argv[1])
info2 = linkerinfoRead(sys.argv[2])
fq1 = FqRead2(sys.argv[3])
fq2 = FqRead2(sys.argv[4])
foutBCOinfo1 = writeFile("{0}.BCO.R1.info".format(sys.argv[5]))
foutBCOinfo2 = writeFile("{0}.BCO.R2.info".format(sys.argv[5]))
foutBCOfq1 = writeFile("{0}.BCO.R1.fastq".format(sys.argv[5]))
foutBCOfq2 = writeFile("{0}.BCO.R2.fastq".format(sys.argv[5]))
foutBCAinfo1 = writeFile("{0}.BCA.R1.info".format(sys.argv[5]))
foutBCAinfo2 = writeFile("{0}.BCA.R2.info".format(sys.argv[5]))
foutBCAfq1 = writeFile("{0}.BCA.R1.fastq".format(sys.argv[5]))
foutBCAfq2 = writeFile("{0}.BCA.R2.fastq".format(sys.argv[5]))
for (flag1, outlist1) , (flag2, outlist2), (rid1,rseq1,rsyb1,rqual1), (rid2,rseq2,rsyb2,rqual2) in zip(info1, info2, fq1, fq2):
    info1 = outlist1.split("\t")
    info2 = outlist2.split("\t")
    flag1 = True
    if info1[1] == "-1" or int(info1[2]) > 5:
        flag1 = False
    flag2 = True
    if info2[1] == "-1" or int(info2[2]) > 5:
        flag2 = False
    if flag1 == False and flag2 == True:
        if info2[7] == "BCO":
            foutBCOinfo1.write("{0}\n".format(outlist1))
            foutBCOinfo2.write("{0}\n".format(outlist2))
            foutBCOfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            foutBCOfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
        else:
            foutBCAinfo1.write("{0}\n".format(outlist1))
            foutBCAinfo2.write("{0}\n".format(outlist2))
            foutBCAfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            foutBCAfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
    elif flag1 == True and flag2 == False:
        if info1[7] == "BCO":
            foutBCOinfo1.write("{0}\n".format(outlist1))
            foutBCOinfo2.write("{0}\n".format(outlist2))
            foutBCOfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            foutBCOfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
        else:
            foutBCAinfo1.write("{0}\n".format(outlist1))
            foutBCAinfo2.write("{0}\n".format(outlist2))
            foutBCAfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            foutBCAfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
    elif flag1 == True and flag2 == True:
        if info1[7] == "BCA" and info2[7] == "BCA":
            foutBCAinfo1.write("{0}\n".format(outlist1))
            foutBCAinfo2.write("{0}\n".format(outlist2))
            foutBCAfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            foutBCAfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
        else:
            foutBCOinfo1.write("{0}\n".format(outlist1))
            foutBCOinfo2.write("{0}\n".format(outlist2))
            foutBCOfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            foutBCOfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
foutBCOinfo1.close()
foutBCOinfo2.close()
foutBCOfq1.close()
foutBCOfq2.close()
foutBCAinfo1.close()
foutBCAinfo2.close()
foutBCAfq1.close()
foutBCAfq2.close()