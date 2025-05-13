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

fin = readFile(sys.argv[1])
n = 0
barcode = {}
for line in fin:
    n += 1
    barcode["{0}".format(n)] = line.strip()
fin.close()

fin = readFile(sys.argv[2])
n = 0
feature = {}
for line in fin:
    n += 1
    feature["{0}".format(n)] = line.split()[0]
fin.close()

fin = readFile(sys.argv[3])
fin.readline()
fin.readline()
fin.readline()
for line in fin:
    tmp = line.strip().split()
    print("{0}\t{1}\t{2}".format(feature[tmp[0]],barcode[tmp[1]],tmp[2]))
fin.close()