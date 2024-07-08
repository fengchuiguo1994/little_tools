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
    parser.add_argument("-r", "--report",      required=True, help="the report file")
    parser.add_argument('-g', "--gemome",      required=True, help="The genome file")
    parser.add_argument('-t', "--gtffile",     required=True, help="The gtf file")
    parser.add_argument('-f', "--gfffile",     required=True, help="The gff file")
    parser.add_argument('-og', "--outgemome",  required=True, help="The output genome file")
    parser.add_argument('-ot', "--outgtffile", required=True, help="The output gtf file")
    parser.add_argument('-of', "--outgfffile", required=True, help="The output gff file")
    args = parser.parse_args()

    mydict = {}
    fin = readFile(args.report)
    for line in fin:
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        if tmp[-1].startswith("chrUn_") or tmp[-1].endswith("_random"):
            continue
        mydict[tmp[6]] = tmp[-1]
    fin.close()
    print(mydict)
    sys.exit(1)

    fout = writeFile(args.outgemome)
    fin = FaRead(args.gemome)
    for fid,fseq in fin:
        if fid in mydict:
            fout.write(">{0}\n".format(mydict[fid]))
            fout.write("{0}\n".format(fseq))
    fout.close()

    fout = writeFile(args.outgtffile)
    fin = readFile(args.gtffile)
    for line in fin:
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        if tmp[2] == "region":
            continue
        if tmp[0] in mydict:
            fout.write("{0}\t{1}\n".format(mydict[tmp[0]], "\t".join(tmp[1:])))
    fout.close()
    
    fout = writeFile(args.outgfffile)
    fin = readFile(args.gfffile)
    for line in fin:
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        if tmp[2] == "region":
            continue
        if tmp[0] in mydict:
            fout.write("{0}\t{1}\n".format(mydict[tmp[0]], "\t".join(tmp[1:])))
    fout.close()
        