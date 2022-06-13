import sys
from math import ceil
import re
import argparse
import subprocess
import numpy as np

parser = argparse.ArgumentParser(description="get the interaction from bedpe file")
parser.add_argument("-p","--bedpe",required=True,help="the input bedpe file")
parser.add_argument("-g","--genome",required=True,help="the genome size file")
parser.add_argument("-o","--outputprefix",default="out",type=str,help="the output prefix [out]")
parser.add_argument("-b","--binsize",default=1000000,type=int,help="the resolution [1000000]")
parser.add_argument("-c","--combine",default=0,type=int,help="0:a[i][j]!=a[j][i] other:a[i][j]=a[j][i]=(a[i][j]+a[j][i]) [0]")
parser.add_argument("-a","--all",default=1,type=int,help="1:calculator whole genome,other calculator in intra-chromosome [1]")
args = parser.parse_args()

def parser_bed(bed,contain,mydic):
    flag = None
    n = 0
    N = 0
    with open(bed,'r') as fbed:
        for line in fbed:
            if line.strip()=="":
                continue
            N += 1
            tmp = line.strip().split("\t")
            tmp[3] = int(tmp[3])
            if flag == None:
                if tmp[0] not in mydic:
                    mydic[tmp[0]] = {}
                mydic[tmp[0]][tmp[3]-n] = tmp[3]
            else:
                if flag == tmp[0]:
                    mydic[tmp[0]][tmp[3]-n] = tmp[3]
                else:
                    n = tmp[3]-1
                    if tmp[0] not in mydic:
                        mydic[tmp[0]] = {}
                        mydic[tmp[0]][tmp[3]-n] = tmp[3]
            flag = tmp[0]
    for i in mydic.keys():
        contain.add(i)
    return N

def parser_bedpe(bedpe,binsize,all,combine,contain,mat):
    with open(bedpe,'r') as fbedpe:
        for line in fbedpe:
            if line.strip()=="":
                continue
            tmp = line.strip().split("\t")
            if all != 1 and tmp[0]!= tmp[3]:
                continue
            if tmp[0] not in contain or tmp[3] not in contain:
                # print(tmp[0],tmp[3])
                continue
            a = ceil((int(tmp[1])+int(tmp[2]))/2/binsize) # before
            a = mydic[tmp[0]][a]
            b = ceil((int(tmp[4])+int(tmp[5]))/2/binsize) # after
            b = mydic[tmp[3]][b]
            mat[a-1][b-1] += 1
            if combine != 0:
                mat[b-1][a-1] += 1



if __name__ == "__main__":
    returncode,returnresult = subprocess.getstatusoutput("bedtools makewindows -g {0} -w {1} -s {2} | awk -v OFS=\"\\t\" \'{{print $0,NR}}' > {3}.tmp.bed".format(args.genome,args.binsize,args.binsize,args.outputprefix))
    if returncode != 0:
        print ("[ERROR]: the bedtools maybe have some error : {0}\n".format(returnresult))
        exit()
    
    mydic = {}
    contain = set()
    N = parser_bed("{0}.tmp.bed".format(args.outputprefix),contain,mydic)

    mat = np.zeros((N,N),dtype=int)
    parser_bedpe(args.bedpe,args.binsize,args.all,args.combine,contain,mat)

    np.savetxt("{0}.mat".format(args.outputprefix),mat,fmt='%d',delimiter='\t')