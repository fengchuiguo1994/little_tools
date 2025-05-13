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

myhash = {}
for index, infile in enumerate(sys.argv[1].split(",")):
    fin = readFile(infile)
    for line in fin:
        tmp = line.strip().split()
        flag = "\t".join(tmp[0:6])
        val = "\t".join(tmp[6:])
        if flag not in myhash:
            myhash[flag] = []
        while len(myhash[flag]) < index:
            myhash[flag].append("-\t-\t-\t-")
        myhash[flag].append(val)
    fin.close()

samples = index + 1
fout = writeFile(sys.argv[2])
for k in myhash.keys():
    while len(myhash[k]) < samples:
        myhash[k].append("-\t-\t-\t-")
    fout.write("{0}\t{1}\n".format(k,"\t".join(myhash[k])))
fout.close()