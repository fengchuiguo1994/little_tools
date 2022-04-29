### parse the string which is the nian col in gtf file  ###
from collections import OrderedDict
def parseString(mystr):
    """
    mystr: gtf attribution information string
    return: a dictionary
    """
    d = OrderedDict()
    if mystr.endswith(";"):
        mystr = mystr[:-1]
    for i in mystr.strip().split(";"):
        tmp = i.strip().split()
        if tmp == "":
            continue
        d[tmp[0]] = tmp[1].replace("\"","")
    return d

### Operating common files  ###
import gzip
def readFile(infile):
    """
    infile: input file
    return: file handle
    """
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
        
def writeFile(outfile):
    """
    outfile: output file
    return: file handle
    """
    if outfile.endswith((".gz","gzip")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

### Operating sam/bam files  ###
import pysam
def readSam(insamfile):
    """
    insamfile: input sam/bam file
    return: file handle
    """
    if insamfile.endswith(".bam"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam.gz"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam"):
        insam = pysam.AlignmentFile(insamfile,'r')
    else:
        raise ValueError("the input sam/bam file is not end with sam or bam!")
    return insam
        
def writeSam(outsamfile,header):
    """
    outsamfile: output sam/bam file
    header: the sam/bam file's header(chromosome information, created by insam.handle)
    return: file handle
    """
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

class SamRead:
    """
    get sam read
    """
    def __init__(self,infile):
        self.fin = readSam(infile)
    def __iter__(self):
        flag = None
        outlist = []
        for read in self.fin:
            if flag != None and flag != read.query_name:
                yield outlist
                outlist = []
            outlist.append(read)
            flag = read.query_name
        self.fin.close()
        yield outlist
    def header(self):
        return self.fin.header

import sys
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
            elif line.startswith(">"):
                yield flag,outstr
                flag = line.strip()[1:]
                outstr = ""
            else:
                outstr += line.strip()
        self.fin.close()
        yield flag,outstr

class FqRead:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip().split()[0]
            rid = rid[1:]
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()

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
        
def rc(sequence):
    """
    Reverse complementary sequence
    sequence: 
    return: 
    """
    seq = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = seq.translate(trantab)
    return string

### Interval tree ###
from bx.intervals.intersection import Intersecter, Interval
def regionTree(tmp,resFrag):
    """
    tmp: The length of the list is at least 3(chrom,start,end)
    resFrag: a dictionary to store the Interval tree
    """
    if tmp[0] not in resFrag:
        resFrag[tmp[0]] = Intersecter()
    resFrag[tmp[0]].add_interval(Interval(int(tmp[1]),int(tmp[2]),tmp[3:]))
    
def regionFind(tree,start,end):
    """
    tree: the intervals tree
    start: the given region's start
    end: the given region's end
    return: a list of overlap region
    """
    return tree.find(start,end)

### reverse strand ###
def reverseStrand(mystr):
    """
    mystr: the strand( + or -)
    return: the opposite strand(- or +)
    """
    if mystr == "+":
        return '-'
    elif mystr == "-":
        return "+"
    else:
        raise ValueError("please check the input strand, must be +/-\n{0}\n".format(mystr))

### determine whether there is overlap between the two intervals ###
def overlap(c1,s1,e1,c2,s2,e2):
    if c1 == c2:
        return e1 > s2 and s1 < e2
    return False

def detOverlap(s1,e1,s2,e2):
    """
    s1,e1,s2,e2: two intervals start and end pos
    return: the True/False whether two intervals has overlap
    """
    if s1 <= s2:
        if e1 > s2:
            return True
        else:
            return False
    else:
        if e2 > s1:
            return True
        else:
            return False
def detOverlapDis(s1,e1,s2,e2):
    """
    s1,e1,s2,e2: two intervals start and end pos
    return: the True/False whether two intervals has overlap and the overlap length or distance
    """
    if s1 <= s2:
        if e1 > s2:
            if e1 < e2:
                return True,e1-s2
            else:
                return True,e2-s2
        else:
            return False,s2-e1+1
    else:
        if e2 > s1:
            if e2 < e1:
                return True,e2-s1
            else:
                return True,e1-s1
        else:
            return False,s1-e2+1
def detOverlapChrom(c1,s1,e1,c2,s2,e2):
    """
    c1,s1,e1,c2,s2,e2: two intervals chrom start and end pos
    return: the True/False whether two intervals has overlap
    """
    if c1 == c2:
        return detOverlap(s1,e1,s2,e2)
    else:
        return False
def detOverlapDisChrom(c1,s1,e1,c2,s2,e2):
    """
    c1,c2,s1,e1,s2,e2: two intervals chrom start and end pos
    return: the True/False whether two intervals has overlap and the overlap length or distance.
    if local on two different chromosome, the distance is -1
    """
    if c1 == c2:
        return detOverlapDis(s1,e1,s2,e2)
    else:
        return False,-1

### merge the overlap interva ( like bedtools merge ) ###
def mergeRegion(regions,order=True):
    """
    regions: A list of locations [[start1,end1],[start2,end2],[start3,end3]...]
    order: is the regions in order
    return: the merge regions
    """
    import copy
    if len(regions) < 2:
        return copy.deepcopy(regions)
    if order != True:
        tmpr = sorted(regions,key=lambda x:(x[0],x[1]))
    else:
        tmpr = copy.deepcopy(regions)
    merge = []
    start = tmpr[0][0]
    end = tmpr[0][1]
    for s,e in tmpr:
        if s > end:
            merge.append([start,end])
            start = s
            end = e
        else:
            if e > end:
                end = e
    merge.append([start,end])
    return merge
### extract blank area ( like get the intron region from exon regions ) ###
def getBlankRegin(regions,order=True,merge=True):
    """
    regions: A list of locations [[start1,end1],[start2,end2],[start3,end3]...]
    order: is the regions in order
    merge: the region has no overlap
    return: the black regions
    """
    if order != True:
        raise ValueError("the input regions must be sorted")
    if merge != True:
        raise ValueError("the input regions has no overlap")
    black = []
    if len(regions) > 1:
        start = regions[0][1]
        for s,e in regions[1:]:
            black.append([start,s])
            start = e
    return black


### interval class and some method ###
class region():
    """
    A region
    start: a region start pos
    end: a region end pos
    """
    def __init__(self,start,end):
        self.start = int(start)
        self.end = int(end)
        if self.start > self.end:
            raise ValueError("start is larger than end")
    def __repr__(self):
        return "class {0}: {1}\t{2}".format(self.__class__.__name__,self.start,self.end)
    def __str__(self):
        return self.__repr__()
    def __eq__(self,other):
        if not isinstance(other,region):
            return False
        return self.start == other.start and self.end == other.end
    def __neq__(self,other):
        if not isinstance(other,region):
            return True
        return not self.__eq__(other)
    def __len__(self):
        return self.end-self.start
    def length(self):
        return self.__len__

class bed3(region):
    """
    A bed3 format region
    chrom: a chromosme
    start: a region start pos
    end: a region end pos
    """
    def __init__(self,chrom,start,end):
        region.__init__(self,start,end)
        self.chrom = chrom
    def __repr__(self):
        return "class {0}: {1}\t{2}\t{3}".format(self.__class__.__name__,self.chrom,self.start,self.end)
    def __str__(self):
        return self.__repr__()
    def __eq__(self,other):
        if not isinstance(other,bed3):
            return False
        return super(bed3, self).__eq__(other) and self.chrom == other.chrom
    def __neq__(self,other):
        if not isinstance(other,bed3):
            return True
        return not self.__eq__(other)

class bed6(bed3):
    """
    A bed6 format region
    chrom: a chromosome
    start: a region start pos
    end: a region end pos
    name: region name
    score: score
    strand: + or - strand
    """
    def __init__(self,chrom,start,end,name=".",score=0,strand="."):
        bed3.__init__(self,chrom,start,end)
        self.name = name
        self.score = score
        self.strand = strand
    def __repr__(self):
        return "class {0}: {1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(self.__class__.__name__,self.chrom,self.start,self.end,self.name,self.score,self.strand)
    def __str__(self):
        return self.__repr__()
    def __eq__(self,other):
        if not isinstance(other,bed6):
            return False
        return super(bed6, self).__eq__(other) and self.name == other.name and self.score == other.score and self.strand == other.strand
    def __neq__(self,other):
        if not isinstance(other,bed6):
            return True
        return not self.__eq__(other)