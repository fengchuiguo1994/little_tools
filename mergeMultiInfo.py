import sys

import gzip
def readFile(infile):
    if infile.endswith((".gz","gzip",".GZ",".GZIP")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
def writeFile(outfile):
    if outfile.endswith((".gz","gzip",".GZ",".GZIP")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

infiles = sys.argv[1].split(",")
outfile = sys.argv[2]
nfile = len(infiles)

interaction = {}
for index,infile in enumerate(infiles):
    fin = readFile(infile)
    for line in fin:
        tmp = line.strip().split()
        flag = tmp[0]
        if flag not in interaction:
            interaction[flag] = []
        while len(interaction[flag]) < index:
            interaction[flag].append("0")
        interaction[flag].append(tmp[1])
    fin.close()

fout = writeFile(outfile)
for k in interaction.keys():
    while len(interaction[k]) < nfile:
        interaction[k].append("0")
    fout.write("{0}\t{1}\n".format(k,"\t".join(interaction[k])))
fout.close