import sys
import collections

import gzip
def readFile(infile):
    """
    infile: input file
    return: file handle
    """
    if infile.endswith((".gz","gzip",".GZ",".GZIP")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin

def writeFile(outfile):
    """
    outfile: output file
    return: file handle
    """
    if outfile.endswith((".gz","gzip",".GZ",".GZIP")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

fin = readFile(sys.argv[1])
fout = writeFile(sys.argv[2])

mydict = collections.OrderedDict()
next(fin)
for line in fin:
    tmp = line.strip().split()
    if len(tmp) != 9:
        sys.exit("file format error!")
    if tmp[0] not in mydict:
        mydict[tmp[0]] = 0
    mydict[tmp[0]] += float(tmp[7])

for k in mydict.keys():
    fout.write("{0}\t{1}\n".format(k,mydict[k]))

fin.close()
fout.close()