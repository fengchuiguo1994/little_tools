import sys
import pysam

fin = pysam.AlignmentFile(sys.stdin)
print(fin.header,end="")
for read in fin:
    rid,barcode = read.query_name.split("-")
    read.query_name = rid
    read.set_tag("BX",barcode)
    print(read.tostring())