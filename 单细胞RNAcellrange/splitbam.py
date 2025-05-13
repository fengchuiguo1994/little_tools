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

fin = readSam(sys.argv[1])
fout1 = writeSam("test0/possorted_genome_bam.bam", header=fin.header)
fout2 = writeSam("test17/possorted_genome_bam.bam", header=fin.header)
fout3 = writeSam("test19/possorted_genome_bam.bam", header=fin.header)
fout4 = writeSam("test25/possorted_genome_bam.bam", header=fin.header)

for read in fin:
    if read.get_tag("xf") == 0:
        fout1.write(read)
    elif read.get_tag("xf") == 17:
        fout2.write(read)
    elif read.get_tag("xf") == 19:
        fout3.write(read)
    elif read.get_tag("xf") == 25:
        fout4.write(read)
    else:
        print(read)
        sys.exit("error")
fin.close()
fout1.close()
fout2.close()
fout3.close()
fout4.close()