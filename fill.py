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
fout = writeFile(sys.argv[2])
fin.readline()
flagc = ""
flags = 0
for line in fin:
    tmp = line.split(maxsplit=2)
    if flagc != "" and flagc != tmp[0]:
        fout.write("{0}\t{1}\t{2}\t{3}".format(tmp[0],0,tmp[1],tmp[2]))
        flagc = tmp[0]
        flags = tmp[1]
    else:
        fout.write("{0}\t{1}\t{2}\t{3}".format(tmp[0],flags,tmp[1],tmp[2]))
        flagc = tmp[0]
        flags = tmp[1]
fin.close()
fout.close()
