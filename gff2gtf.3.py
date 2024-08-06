import sys
import re

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

def deal(mylist):
    if len(mylist) == 2:
        gid = "MT-{0}".format(re.search(r'description=(.+?);',mylist[0][8])[1])
        print('ChrMt\t{0}\tgene\t{1}\t{2}\t.\t{3}\t.\tgene_id "{4}"; gene_biotype "protein_coding";'.format(mylist[0][1],mylist[0][3],mylist[0][4],mylist[0][6],gid))
        print('ChrMt\t{0}\ttranscript\t{1}\t{2}\t.\t{3}\t.\tgene_id "{4}"; transcript_id "{4}-RNA";'.format(mylist[0][1],mylist[0][3],mylist[0][4],mylist[0][6],gid))
        print('ChrMt\t{0}\texon\t{1}\t{2}\t.\t{3}\t.\tgene_id "{4}"; transcript_id "{4}-RNA";'.format(mylist[0][1],mylist[0][3],mylist[0][4],mylist[0][6],gid))
    elif len(mylist) == 3:
        gid = "MT-{0}".format(re.search(r'ID=(.+?);',mylist[0][8])[1])
        print('ChrMt\t{0}\tgene\t{1}\t{2}\t.\t{3}\t.\tgene_id "{4}"; gene_biotype "{5}";'.format(mylist[0][1],mylist[0][3],mylist[0][4],mylist[0][6],gid,mylist[1][2]))
        print('ChrMt\t{0}\ttranscript\t{1}\t{2}\t.\t{3}\t.\tgene_id "{4}"; transcript_id "{4}-RNA";'.format(mylist[0][1],mylist[0][3],mylist[0][4],mylist[0][6],gid))
        print('ChrMt\t{0}\texon\t{1}\t{2}\t.\t{3}\t.\tgene_id "{4}"; transcript_id "{4}-RNA";'.format(mylist[0][1],mylist[0][3],mylist[0][4],mylist[0][6],gid))
    
fin = readFile(sys.argv[1])
flag = None
mylist = []
for line in fin:
    if line.startswith("#"):
        continue
    tmp = line.strip().split("\t")
    if tmp[2] == "region" or tmp[2] == "origin_of_replication" or tmp[2] == "D_loop":
        continue
    if flag != None and tmp[2] == "gene":
        deal(mylist)
        mylist = []
    flag = tmp[2]
    mylist.append(tmp)
deal(mylist)