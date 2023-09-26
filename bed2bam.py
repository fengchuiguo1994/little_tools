import sys
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
fout = pysam.AlignmentFile(sys.argv[3], "wb", template=fin)
fin.close()
fin = readFile(sys.argv[2])
index = 0
for read in fin:
    index += 1
    tmp = read.strip().split()
    start = int(tmp[1])
    end = int(tmp[2])
    a = pysam.AlignedSegment()
    a.query_name = "read{0}".format(index)
    a.query_sequence=""
    a.query_qualities = pysam.qualitystring_to_array("")
    a.flag = 8+32+64
    a.reference_id = fout.get_tid(tmp[0])
    a.reference_start = start
    a.mapping_quality = 60
    a.cigar = ((0,end-start),)
    # AS:i:81 XS:i:19 BX:Z:AGACAGGCAGACTGTT
    a.tags = (("AS", end-start), ("XS", 0),("BX", tmp[3]),("CB", tmp[3]))
    fout.write(a)
fout.close()