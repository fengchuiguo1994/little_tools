import sys

from collections import OrderedDict
def parseString(mystr):
    d = OrderedDict()
    if mystr.endswith(";"):
        mystr = mystr[:-1]
    for i in mystr.strip().split(";"):
        tmp = i.strip().split()
        if tmp == "":
            continue
        d[tmp[0]] = tmp[1].replace("\"","")
    return d

from bx.intervals.intersection import Intersecter, Interval
def regionTree(tmp,resFrag):
    if tmp[0] in resFrag:
        resFrag[tmp[0]].add_interval(Interval(int(tmp[1]),int(tmp[2])))
    else:
        resFrag[tmp[0]] = Intersecter()
        resFrag[tmp[0]].add_interval(Interval(int(tmp[1]),int(tmp[2])))
def regionTree1(mylist):
    resFrag = Intersecter()
    for s,e in mylist:
        resFrag.add_interval(Interval(s,e))
    return resFrag
def regionTree2(tmp,resFrag,strand,info,geneid):
    if strand not in resFrag:
        resFrag[strand] = {}
    if tmp[0] in resFrag[strand]:
        resFrag[strand][tmp[0]].add_interval(Interval(int(tmp[1]),int(tmp[2]),value={"exon":info,'geneid':geneid}))
    else:
        resFrag[strand][tmp[0]] = Intersecter()
        resFrag[strand][tmp[0]].add_interval(Interval(int(tmp[1]),int(tmp[2]),value={"exon":info,'geneid':geneid}))
def regionFind(tree,start,end):
    return tree.find(start,end)

import pysam
def openSam(insamfile):
    if insamfile.endswith(".bam"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam"):
        insam = pysam.AlignmentFile(insamfile,'r')
    else:
        raise ValueError("the input sam/bam file is not end with sam or bam!")
    return insam
def writeSam(outsamfile,header):
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

def reverseStrand(mystr):
    if mystr == "+":
        return '-'
    elif mystr == "-":
        return "+"
    else:
        raise ValueError("please check the input strand, must be +/-\n{0}\n".format(mystr))

import gzip
def openFile(infile):
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

class GeneTrack:
    def __init__(self,geneid,genename=None,chrom=None,strand=None,exon=[]):
        self.geneid = geneid
        self.genename = genename
        self.chrom = chrom
        self.strand = strand
        self.exon = exon

    def getMergeExon(self):
        self.exon.sort(key=lambda x:(x[0],x[1]))
        self.mergeExon = []
        start = end = -1
        for s,e in self.exon:
            if start == -1:
                start = s
                end = e
            else:
                if s > end:
                    self.mergeExon.append([start,end])
                    start = s
                    end = e
                else:
                    if e > end:
                        end = e
        self.mergeExon.append([start,end])
        self.genestart = self.mergeExon[0][0]
        self.geneend = self.mergeExon[-1][-1]

def getGeneTree(ingtf):
    fin = openFile(ingtf)
    resFrag = {} # return
    GUcount = {} # return
    GLcount = {} # return
    EUcount = {} # return
    ELcount = {} # return
    geneinfo = {}
    for line in fin:
        if isinstance(line,bytes):
            line = line.decode()
        if line == "\n":
            continue
        if line.startswith("#"):
            continue
        (seqname,source,feature,start,end,score,strand,frame,attribution) = line.strip().split("\t",8)
        if feature == "exon":
            attr = parseString(attribution)
            if "gene_id" not in attr:
                sys.exit("some record has no gene_id info, please check your gtf file\n")
            if strand == ".":
                strand = "+"
                sys.stderr.write("{0}'s strand is \'.\' not be \'+\' or \'-\'. we treat it to \'+\'".format(attr["gene_id"]))
            if "gene_name" not in attr:
                attr["gene_name"] = attr["gene_id"]
            if attr["gene_id"] not in geneinfo:
                geneinfo[attr["gene_id"]] = GeneTrack(geneid=attr["gene_id"],genename=attr["gene_name"],chrom=seqname,strand=strand,exon=[])
            geneinfo[attr["gene_id"]].exon.append([int(start)-1,int(end)])
            if strand != geneinfo[attr["gene_id"]].strand:
                raise ValueError("please check the gtf file, a gene has two strand info\n{0}\n".format(line))
            if seqname != geneinfo[attr["gene_id"]].chrom:
                raise ValueError("please check the gtf file, a gene has two chromosome info\n{0}\n".format(line))
    for i in geneinfo.keys():
        geneinfo[i].getMergeExon()
        regionTree2([geneinfo[i].chrom,geneinfo[i].genestart,geneinfo[i].geneend],resFrag,geneinfo[i].strand,regionTree1(geneinfo[i].mergeExon),geneinfo[i].geneid)
        GUcount[geneinfo[i].geneid] = 0
        GLcount[geneinfo[i].geneid] = 0
        EUcount[geneinfo[i].geneid] = 0
        ELcount[geneinfo[i].geneid] = 0
    fin.close()
    return resFrag,GUcount,GLcount,EUcount,ELcount

def detOverlapDis(s1,e1,s2,e2):
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

def countlength(overlapgenelengthdict,start,end,mylist):
    for i in mylist:
        if i.value['geneid'] not in overlapgenelengthdict:
            overlapgenelengthdict[i.value['geneid']] = 0
        ovf, dis = detOverlapDis(start,end,i.start,i.end)
        if not ovf:
            sys.exit("may have some bug1, interval tree overlap has some bug")
        overlapgenelengthdict[i.value['geneid']] += dis

def getOverlapGene(read,flag,genetree,strand,chrom):
    overlapgene = []
    overlapgenelengthdict = {}
    for start,end in read.blocks:
        if flag == "yes":
            tmp = regionFind(genetree[strand][chrom],start,end)
            overlapgene.extend(tmp)
            countlength(overlapgenelengthdict,start,end,tmp)
        elif flag == "reverse":
            rcstrand = reverseStrand(strand)
            tmp = regionFind(genetree[rcstrand][chrom],start,end)
            overlapgene.extend(tmp)
            countlength(overlapgenelengthdict,start,end,tmp)
        else:
            tmp = regionFind(genetree[strand][chrom],start,end)
            overlapgene.extend(tmp)
            countlength(overlapgenelengthdict,start,end,tmp)
            rcstrand = reverseStrand(strand)
            tmp = regionFind(genetree[rcstrand][chrom],start,end)
            overlapgene.extend(tmp)
            countlength(overlapgenelengthdict,start,end,tmp)
    overlapgeneid = set()
    tmpoverlapgene = []
    for i in overlapgene:
        if i.value['geneid'] not in overlapgeneid:
            overlapgeneid.add(i.value['geneid'])
            tmpoverlapgene.append(i)
    return tmpoverlapgene,overlapgenelengthdict

def getOverlapExon(read,flag,tree):
    overlapexon = []
    length = 0
    for start,end in read.blocks:
        tmp = regionFind(tree,start,end)
        overlapexon.extend(tmp)
        for i in tmp:
            ovf, dis = detOverlapDis(start,end,i.start,i.end)
            if not ovf:
                sys.exit("may have some bug2, interval tree overlap has some bug")
            length += dis
    return overlapexon,length

def topGene(overlapgenelengthdict):
    mytmp = sorted(overlapgenelengthdict.items(),key=lambda x:x[1],reverse=True)
    aa = mytmp[0][1]
    outlist = []
    for i,j in mytmp:
        if j > aa-15:
            outlist.append(i)
        else:
            break
    return outlist

def main(samfile,ingtf,outsamfile,flag):
    markstrand = ("GU","GL","EU","EL")
    marktype = ("Unmapped","NoFeatures","Features","Ambiguity")
    genetree,GUcount,GLcount,EUcount,ELcount = getGeneTree(ingtf)
    othercount = {}
    for i in markstrand:
        othercount[i] = {}
        for j in marktype:
            othercount[i][j] = 0
    insam = openSam(samfile)
    outsam = writeSam(outsamfile,insam.header)

    for read in insam:
        if read.is_unmapped:
            read.set_tag("GU",'Unmapped')
            read.set_tag("GL",'Unmapped')
            read.set_tag("EU",'Unmapped')
            read.set_tag("EL",'Unmapped')
            othercount["GU"]["Unmapped"] += 1
            othercount["GL"]["Unmapped"] += 1
            othercount["EU"]["Unmapped"] += 1
            othercount["EL"]["Unmapped"] += 1
            outsam.write(read)
            continue
        chrom = read.reference_name
        if read.is_reverse:
            strand = "-"
        else:
            strand = "+"

        ### overlap gene region
        overlapgene,overlapgenelengthdict = getOverlapGene(read,flag,genetree,strand,chrom)

        ### overlap exon region
        if len(overlapgene) == 0:
            read.set_tag("GU",'NoFeatures')
            read.set_tag("GL",'NoFeatures')
            read.set_tag("EU",'NoFeatures')
            read.set_tag("EL",'NoFeatures')
            othercount["GU"]["NoFeatures"] += 1
            othercount["GL"]["NoFeatures"] += 1
            othercount["EU"]["NoFeatures"] += 1
            othercount["EL"]["NoFeatures"] += 1
            outsam.write(read)
        elif len(overlapgene) == 1:
            geneid = overlapgene[0].value['geneid']
            read.set_tag("GU",geneid)
            read.set_tag("GL",geneid)
            GUcount[geneid] += 1
            GLcount[geneid] += 1
            othercount["GU"]["Features"] += 1
            othercount["GL"]["Features"] += 1
            overlapexon,exonlength = getOverlapExon(read,flag,overlapgene[0].value['exon'])
            if len(overlapexon) > 0:
                read.set_tag("EU",geneid)
                read.set_tag("EL",geneid)
                EUcount[geneid] += 1
                ELcount[geneid] += 1
                othercount["EU"]["Features"] += 1
                othercount["EL"]["Features"] += 1
            else:
                read.set_tag("EU",'NoFeatures')
                read.set_tag("EL",'NoFeatures')
                othercount["EU"]["NoFeatures"] += 1
                othercount["EL"]["NoFeatures"] += 1
            outsam.write(read)
        else:
            geneid = [x.value["geneid"] for x in overlapgene]
            read.set_tag("GU", "Ambiguity:[{0}]".format(",".join(geneid)))
            othercount["GU"]["Ambiguity"] += 1
            geneid = topGene(overlapgenelengthdict)
            if len(geneid) > 1:
                read.set_tag("GL", "Ambiguity:[{0}]".format(",".join(geneid)))
                othercount["GL"]["Ambiguity"] += 1
            elif len(geneid) == 1:
                read.set_tag("GL", geneid[0])
                othercount["GL"]["Features"] += 1
                GLcount[geneid[0]] += 1
            else:
                sys.exit("have some bug1, 0 gene has mark, Please feedback to the author\n")

            overlapexonlengthdict = {}
            for ii in overlapgene:
                overlapexon,exonlength = getOverlapExon(read,flag,ii.value['exon'])
                if len(overlapexon) > 0:
                    overlapexonlengthdict[ii.value["geneid"]] = exonlength
            if len(overlapexonlengthdict) == 0:
                read.set_tag("EU",'NoFeatures')
                read.set_tag("EL",'NoFeatures')
                othercount["EU"]["NoFeatures"] += 1
                othercount["EL"]["NoFeatures"] += 1
            elif len(overlapexonlengthdict) == 1:
                exonid = list(overlapexonlengthdict.keys())
                read.set_tag("EU",exonid[0])
                read.set_tag("EL",exonid[0])
                othercount["EU"]["Features"] += 1
                othercount["EL"]["Features"] += 1
                EUcount[exonid[0]] += 1
                ELcount[exonid[0]] += 1
            else:
                exonid = list(overlapexonlengthdict.keys())
                read.set_tag("EU","Ambiguity:[{0}]".format(",".join(exonid)))
                othercount["EU"]["Ambiguity"] += 1
                exonid = topGene(overlapexonlengthdict)
                if len(exonid) > 1:
                    read.set_tag("EL", "Ambiguity:[{0}]".format(",".join(geneid)))
                    othercount["EL"]["Ambiguity"] += 1
                elif len(exonid) == 1:
                    read.set_tag("EL", exonid[0])
                    othercount["EL"]["Features"] += 1
                    ELcount[exonid[0]] += 1
                else:
                    sys.exit("have some bug2, 0 gene has mark, Please feedback to the author\n")
            outsam.write(read)

    insam.close()
    outsam.close()
    fout = writeFile("{0}.raw.count.txt".format(outsamfile))
    for i in sorted(GUcount.keys()):
        fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(i,GUcount[i],GLcount[i],EUcount[i],ELcount[i]))
    fout.close()
    fout = writeFile("{0}.log".format(outsamfile))
    for i in marktype:
        fout.write(i)
        for j in markstrand:
            fout.write("\t{0}".format(othercount[j][i]))
        fout.write("\n")
    fout.close()

if __name__ == "__main__":
    # sam/bam file, gtf file, outfile, flag
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
