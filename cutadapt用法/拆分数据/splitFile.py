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
                yield flag,outlist
                flag = line.strip().split()[0]
                outlist = []
            outlist.append(line.strip())
        self.fin.close()
        yield flag,outlist

info1 = linkerinfoRead(sys.argv[1])
info2 = linkerinfoRead(sys.argv[2])
fq1 = FqRead2(sys.argv[3])
fq2 = FqRead2(sys.argv[4])
fout22info1 = writeFile("{0}.22.R1.info".format(sys.argv[5]))
fout22info2 = writeFile("{0}.22.R2.info".format(sys.argv[5]))
fout22fq1 = writeFile("{0}.22.R1.fastq".format(sys.argv[5]))
fout22fq2 = writeFile("{0}.22.R2.fastq".format(sys.argv[5]))
fout21info1 = writeFile("{0}.21.R1.info".format(sys.argv[5]))
fout21info2 = writeFile("{0}.21.R2.info".format(sys.argv[5]))
fout21fq1 = writeFile("{0}.21.R1.fastq".format(sys.argv[5]))
fout21fq2 = writeFile("{0}.21.R2.fastq".format(sys.argv[5]))
fout12info1 = writeFile("{0}.12.R1.info".format(sys.argv[5]))
fout12info2 = writeFile("{0}.12.R2.info".format(sys.argv[5]))
fout12fq1 = writeFile("{0}.12.R1.fastq".format(sys.argv[5]))
fout12fq2 = writeFile("{0}.12.R2.fastq".format(sys.argv[5]))
fout11dropinfo1 = writeFile("{0}.11.drop.R1.info".format(sys.argv[5]))
fout11dropinfo2 = writeFile("{0}.11.drop.R2.info".format(sys.argv[5]))
fout11dropfq1 = writeFile("{0}.11.drop.R1.fastq".format(sys.argv[5]))
fout11dropfq2 = writeFile("{0}.11.drop.R2.fastq".format(sys.argv[5]))
fout11shiftinfo1 = writeFile("{0}.11.shift.R1.info".format(sys.argv[5]))
fout11shiftinfo2 = writeFile("{0}.11.shift.R2.info".format(sys.argv[5]))
fout11shiftfq1 = writeFile("{0}.11.shift.R1.fastq".format(sys.argv[5]))
fout11shiftfq2 = writeFile("{0}.11.shift.R2.fastq".format(sys.argv[5]))
fout11fltinfo1 = writeFile("{0}.11.flt.R1.info".format(sys.argv[5]))
fout11fltinfo2 = writeFile("{0}.11.flt.R2.info".format(sys.argv[5]))
fout11fltfq1 = writeFile("{0}.11.flt.R1.fastq".format(sys.argv[5]))
fout11fltfq2 = writeFile("{0}.11.flt.R2.fastq".format(sys.argv[5]))
for (flag1, outlist1) , (flag2, outlist2), (rid1,rseq1,rsyb1,rqual1), (rid2,rseq2,rsyb2,rqual2) in zip(info1, info2, fq1, fq2):
    if len(outlist1) > 1 and len(outlist2) > 1:
        for k in outlist1:
            fout22info1.write("{0}\n".format(k))
        for k in outlist2:
            fout22info2.write("{0}\n".format(k))
        fout22fq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
        fout22fq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
    elif len(outlist1) > 1:
        for k in outlist1:
            fout21info1.write("{0}\n".format(k))
        for k in outlist2:
            fout21info2.write("{0}\n".format(k))
        fout21fq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
        fout21fq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
    elif len(outlist2) > 1:
        for k in outlist1:
            fout12info1.write("{0}\n".format(k))
        for k in outlist2:
            fout12info2.write("{0}\n".format(k))
        fout12fq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
        fout12fq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
    else:
        info1 = outlist1[0].split("\t")
        info2 = outlist2[0].split("\t")
        if info1[1] == "-1" or info2[1] == "-1":
            for k in outlist1:
                fout11dropinfo1.write("{0}\n".format(k))
            for k in outlist2:
                fout11dropinfo2.write("{0}\n".format(k))
            fout11dropfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            fout11dropfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
            continue
        flag1 = True
        if int(info1[2]) > 5:
            flag1 = False
        flag1type = info1[7]
        flag2 = True
        if int(info2[2]) > 5:
            flag2 = False
        flag2type = info2[7]
        if flag1 == False or flag2 == False:
            for k in outlist1:
                fout11shiftinfo1.write("{0}\n".format(k))
            for k in outlist2:
                fout11shiftinfo2.write("{0}\n".format(k))
            fout11shiftfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            fout11shiftfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
        else:
            for k in outlist1:
                fout11fltinfo1.write("{0}\n".format(k))
            for k in outlist2:
                fout11fltinfo2.write("{0}\n".format(k))
            fout11fltfq1.write("{0}\n{1}\n{2}\n{3}\n".format(rid1,rseq1,rsyb1,rqual1))
            fout11fltfq2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))