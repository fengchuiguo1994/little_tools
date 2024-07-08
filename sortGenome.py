import argparse
import sys
import random
import copy

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

class FaRead:
    def __init__(self,infile):
        if infile == None or infile == "-":
            self.fin = sys.stdin
        else:
            self.fin = readFile(infile)
    def __iter__(self):
        flag = None
        outstr = ""
        for line in self.fin:
            if flag == None:
                flag = line.strip()[1:]
                flag = flag.split()[0]
            elif line.startswith(">"):
                yield flag,outstr
                flag = line.strip()[1:]
                flag = flag.split()[0]
                outstr = ""
            else:
                outstr += line
        self.fin.close()
        yield flag,outstr

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', "--gemome",      required=True, help="The genome file")
    parser.add_argument('-og', "--outgemome",  required=True, help="The output genome file")
    args = parser.parse_args()

    mydict = {}
    fout = writeFile(args.outgemome)
    fin = FaRead(args.gemome)
    for fid,fseq in fin:
        mydict[fid] = fseq
    for k in sorted(mydict.keys()):
        fout.write(">{0}\n".format(k))
        fout.write("{0}\n".format(mydict[k]))
    fout.close()