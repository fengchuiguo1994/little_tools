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
myhash = {}
for line in fin:
    tmp = line.strip().split()
    if tmp[0] not in myhash:
        myhash[tmp[0]] = writeFile("{0}.{1}.bed".format(sys.argv[2],tmp[0]))
    myhash[tmp[0]].write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(tmp[0], tmp[1], tmp[2], tmp[6], tmp[7], tmp[8]))

    if tmp[3] not in myhash:
        myhash[tmp[3]] = writeFile("{0}.{1}.bed".format(sys.argv[2],tmp[3]))
    myhash[tmp[3]].write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[9]))
fin.close()

for chrom in myhash.keys():
    myhash[chrom].close()
