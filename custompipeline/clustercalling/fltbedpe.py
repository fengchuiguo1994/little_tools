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

fbedpe = readFile(sys.argv[1])
fout = writeFile(sys.argv[2])
for line in fbedpe:
    tmp = line.strip().split()
    if tmp[0] != tmp[3]:
        fout.write(line)
        if tmp[0] > tmp[3]:
            print("1\t{0}".format(line.strip()))
        continue
    if int(tmp[1]) > int(tmp[4]):
        print("2\t{0}".format(line.strip()))
    if int(tmp[5]) - int(tmp[1]) > 8000:
        fout.write(line)
fbedpe.close()
fout.close()