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

class FqRead2:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip()
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--seed", required=False,default=None,type=int, help="the random seed")
    parser.add_argument('-r1', "--read1", required=True, help="The read1 input file")
    parser.add_argument('-r2', "--read2", required=False, default=None, help="The read2 input file")
    parser.add_argument('-o1', "--outputread1", required=True, help="The read1 output file")
    parser.add_argument('-o2', "--outputread2", required=False, default=None, help="The read2 output file, if set, must set the --read2")
    args = parser.parse_args()

    if args.outputread2 != None and args.read2 == None:
        sys.exit("if you give the --outputread2, you must give the read2")

    if args.seed == None:
        num = random.randint(0,99999999999)
        random.seed(num)
        print("current random seed is: {0}".format(num))
    else:
        random.seed(args.seed)

    instrument = ""
    instrument += random.choice("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    instrument += "{0}".format(random.randint(10000,99999))
    runnumber = ""
    runnumber += "{0}".format(random.randint(100,999))
    flowcell = ""
    flowcell += "".join(random.choices("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789",k=9))
    lane = "{0}".format(random.randint(1,9))

    tmp = []
    for i in range(1000,10000):
        tmp.append("{0}".format(i))
    random.shuffle(tmp)
    title = copy.deepcopy(tmp[:1000])

    tmp = []
    for i in range(1000,100000):
        tmp.append("{0}".format(i))
    random.shuffle(tmp)
    xpos = copy.deepcopy(tmp[:4000])

    tmp = []
    for i in range(1000,100000):
        tmp.append("{0}".format(i))
    random.shuffle(tmp)
    ypos = copy.deepcopy(tmp[:4000])

    if args.outputread2 != None: # paired end reads
        fin1 = FqRead2(args.read1)
        fin2 = FqRead2(args.read2)
        fout1 = writeFile(args.outputread1)
        fout2 = writeFile(args.outputread2)
        titleindex = 0
        xposindex = 0
        yposindex = 0
        for index,(r1,r2) in enumerate(zip(fin1,fin2)):
            yposindex += 1
            if yposindex >= 4000:
                yposindex -= 4000
                xposindex += 1
            if xposindex >= 4000:
                xposindex -= 4000
                titleindex += 1
            fout1.write("@{0}:{1}:{2}:{3}:{4}:{5}:{6} 1\n".format(
                instrument,runnumber,flowcell,lane,title[titleindex],xpos[xposindex],ypos[yposindex]))
            fout1.write("{0}\n".format(r1[1]))
            fout1.write("{0}\n".format(r1[2]))
            fout1.write("{0}\n".format(r1[3]))
            fout2.write("@{0}:{1}:{2}:{3}:{4}:{5}:{6} 2\n".format(
                instrument,runnumber,flowcell,lane,title[titleindex],xpos[xposindex],ypos[yposindex]))
            fout2.write("{0}\n".format(r2[1]))
            fout2.write("{0}\n".format(r2[2]))
            fout2.write("{0}\n".format(r2[3]))
        fout1.close()
        fout2.close()
    else:
        fin1 = FqRead2(args.read1)
        fout1 = writeFile(args.outputread1)
        titleindex = 0
        xposindex = 0
        yposindex = 0
        for index,r1 in enumerate(fin1):
            yposindex += 1
            if yposindex >= 4000:
                yposindex -= 4000
                xposindex += 1
            if xposindex >= 4000:
                xposindex -= 4000
                titleindex += 1
            fout1.write("@{0}:{1}:{2}:{3}:{4}:{5}:{6}\n".format(
                instrument,runnumber,flowcell,lane,title[titleindex],xpos[xposindex],ypos[yposindex]))
            fout1.write("{0}\n".format(r1[1]))
            fout1.write("{0}\n".format(r1[2]))
            fout1.write("{0}\n".format(r1[3]))
        fout1.close()