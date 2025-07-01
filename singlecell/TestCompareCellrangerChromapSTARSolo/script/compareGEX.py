import sys
import collections

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

mydict = collections.OrderedDict()
fin1 = readSam(sys.argv[1]) # starsolo
fin2 = readSam(sys.argv[2]) # cellranger
fout = writeFile(sys.argv[3]) # output
for read in fin1:
    if read.is_unmapped:
        flag = "unmapped"
    else:
        flag = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tstarsolo".format(read.query_name, read.reference_name, read.reference_start, read.reference_end, read.mapping_quality, read.get_tag("CB"), read.get_tag("UB"))
    if read.query_name not in mydict:
        mydict[read.query_name] = []
    mydict[read.query_name].append(flag)
fin1.close()

for read in fin2:
    if read.is_unmapped:
        flag = "unmapped"
    else:
        if read.has_tag("CB"):
            cb = read.get_tag("CB")
        else:
            cb = "-"
        if read.has_tag("UB"):
            ub = read.get_tag("UB")
        else:
            ub = "-"
        flag = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tcellranger".format(read.query_name, read.reference_name, read.reference_start, read.reference_end, read.mapping_quality, cb, ub)
    if read.query_name not in mydict:
        mydict[read.query_name] = []
    mydict[read.query_name].append(flag)
fin2.close()

for k in mydict.keys():
    fout.write("{0}\n".format("\n".join(mydict[k])))
fout.close()