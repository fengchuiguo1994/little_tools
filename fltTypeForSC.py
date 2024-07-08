import sys
import gzip
import re

def readFile(infile):
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin

class GTFRead:
    def __init__(self,infile):
        if infile == None or infile == "-":
            self.fin = sys.stdin
        else:
            self.fin = readFile(infile)
    def __iter__(self):
        flag = None
        outstr = []
        for line in self.fin:
            if line.startswith("#"):
                continue
            tmp = line.strip().split("\t")
            if flag == None:
                flag = 2
            elif tmp[2] == "gene":
                yield outstr
                flag = 2
                outstr = []
            outstr.append(line)
        self.fin.close()
        yield outstr

myset = set(["protein_coding", "lncRNA", "C_region", "other", "V_segment"])
fin = GTFRead(sys.argv[1])
for aa in fin:
    # print(aa)
    geneid = re.search(r'''gene_id "(.+?)";''', aa[0])[1]
    # print(geneid)
    for ii in aa:
        if geneid != re.search(r'''gene_id "(.+?)";''', ii)[1]:
            # print(aa)
            sys.exit("check")
    genetype = re.search(r'''gene_biotype "(.+?)";''', aa[0])[1]
    if genetype in myset:
        for ii in aa:
            print(ii, end="")
sys.stdout.flush()
sys.stdout.close()