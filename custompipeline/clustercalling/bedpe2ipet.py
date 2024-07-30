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
    if tmp[0] == tmp[3] and int(tmp[5]) - int(tmp[1]) <= 8000:
        continue
    if tmp[8] == "+":
        start1 = (int(tmp[1]) + int(tmp[2]))//2 - 20
        end1 = start1 + 500
    else:
        end1 = (int(tmp[1]) + int(tmp[2]))//2 + 20
        start1 = end1 - 500

    if tmp[9] == "+":
        start2 = (int(tmp[4]) + int(tmp[5]))//2 - 20
        end2 = start2 + 500
    else:
        end2 = (int(tmp[4]) + int(tmp[5]))//2 + 20
        start2 = end2 - 500
    if tmp[0] == tmp[3] and end2 - start1 <= 8000:
        continue
    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(tmp[0],start1,end1,tmp[3],start2,end2,tmp[6]))
fbedpe.close()
fout.close()