import sys
from bx.intervals.intersection import Intersecter, Interval

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

from bx.intervals.intersection import Intersecter, Interval
def regionTree(tmp,resFrag):
    if tmp[0] not in resFrag:
        resFrag[tmp[0]] = Intersecter()
    resFrag[tmp[0]].add_interval(Interval(tmp[1],tmp[2]))

fanchor = readFile(sys.argv[1])
floop = writeFile(sys.argv[3])
resFrag = {}
for line in fanchor:
    tmp = line.strip().split()
    chrom = tmp[0]
    start = int(tmp[1])
    end = int(tmp[2])
    regionTree([chrom,start,end],resFrag)
fanchor.close()

fbedpe = readFile(sys.argv[2])
covdict = {}
loopdict = {}
for line in fbedpe:
    tmp = line.strip().split()
    r1 = []
    r2 = []
    if tmp[0] in resFrag:
        hit1 = resFrag[tmp[0]].find(int(tmp[1]),int(tmp[2]))
        for k in hit1:
            r1.append("{0}\t{1}\t{2}".format(tmp[0],k.start,k.end))
    if tmp[3] in resFrag:
        hit2 = resFrag[tmp[3]].find(int(tmp[4]),int(tmp[5]))
        for k in hit2:
            r2.append("{0}\t{1}\t{2}".format(tmp[3],k.start,k.end))
    for i in r1:
        for j in r2:
            if i not in loopdict:
                loopdict[i] = {}
            if j not in loopdict[i]:
                loopdict[i][j] = 0
            loopdict[i][j] += 1
    for i in r1:
        if i not in covdict:
            covdict[i] = 0
        covdict[i] += 1
    for i in r2:
        if i not in covdict:
            covdict[i] = 0
        covdict[i] += 1
fbedpe.close()

for k1 in loopdict.keys():
    for k2 in loopdict[k1].keys():
        floop.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k1,k2,loopdict[k1][k2],covdict[k1],covdict[k2],loopdict[k1][k2]/((covdict[k1]*covdict[k2])**0.5)))
floop.close()
