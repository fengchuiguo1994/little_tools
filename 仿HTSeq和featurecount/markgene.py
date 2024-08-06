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

def getGeneTree(ingtf):
    fin = openFile(ingtf)
    resFrag = {} # return
    geneinfo = {}
    genepos = {}
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
            if attr["gene_id"] not in geneinfo:
                geneinfo[attr["gene_id"]] = {}
                genepos[attr["gene_id"]] = {}
                genepos[attr["gene_id"]]['start'] = int(start) - 1
                genepos[attr["gene_id"]]['end'] = int(end)
                genepos[attr["gene_id"]]['strand'] = strand
                genepos[attr["gene_id"]]['chrom'] = seqname
            regionTree([seqname,int(start)-1,int(end)],geneinfo[attr["gene_id"]])
            if int(start) - 1 < genepos[attr["gene_id"]]['start']:
                genepos[attr["gene_id"]]['start'] = int(start) - 1
            if int(end) > genepos[attr["gene_id"]]['end']:
                genepos[attr["gene_id"]]['end'] = int(end)
            if strand != genepos[attr["gene_id"]]['strand']:
                raise ValueError("please check the gtf file, a gene has two strand info\n{0}\n".format(line))
            if seqname != genepos[attr["gene_id"]]['chrom']:
                raise ValueError("please check the gtf file, a gene has two chromosome info\n{0}\n".format(line))
    for i in genepos.keys():
        regionTree2([genepos[i]['chrom'],genepos[i]['start'],genepos[i]['end']],resFrag,genepos[i]['strand'],geneinfo[i],i)
    fin.close()
    return resFrag

def getOverlapGene(read,flag,genetree,strand,chrom):
    overlapgene = []
    for start,end in read.blocks:
        if flag == "yes":
            overlapgene.extend(regionFind(genetree[strand][chrom],start,end))
        elif flag == "reverse":
            rcstrand = reverseStrand(strand)
            overlapgene.extend(regionFind(genetree[rcstrand][chrom],start,end))
        else:
            overlapgene.extend(regionFind(genetree[strand][chrom],start,end))
            rcstrand = reverseStrand(strand)
            overlapgene.extend(regionFind(genetree[rcstrand][chrom],start,end))
    overlapgeneid = set()
    tmpoverlapgene = []
    for i in overlapgene:
        if i.value['geneid'] not in overlapgeneid:
            overlapgeneid.add(i.value['geneid'])
            tmpoverlapgene.append(i)
    return tmpoverlapgene

def getOverlapExon(read,flag,tree):
    overlapexon = []
    for start,end in read.blocks:
        overlapexon.extend(regionFind(tree,start,end))
    return overlapexon

import sys
def main(samfile,ingtf,outsamfile,flag):
    genetree = getGeneTree(ingtf)
    insam = openSam(samfile)
    outsam = writeSam(outsamfile,insam.header)

    for read in insam:
        if read.is_unmapped:
            read.set_tag("SG",'Unmapped')
            read.set_tag("SE",'Unmapped')
            outsam.write(read)
            continue
        chrom = read.reference_name
        if read.is_reverse:
            strand = "-"
        else:
            strand = "+"

        ### overlap gene region
        overlapgene = getOverlapGene(read,flag,genetree,strand,chrom)

        ### overlap exon region
        if len(overlapgene) == 0:
            read.set_tag("SG",'NoFeatures')
            read.set_tag("SE",'NoFeatures')
            outsam.write(read)
        elif len(overlapgene) == 1:
            geneid = overlapgene[0].value['geneid']
            read.set_tag("SG",geneid)
            overlapexon = getOverlapExon(read,flag,overlapgene[0].value['exon'][chrom])
            if len(overlapexon) > 0:
                read.set_tag("SE",geneid)
            else:
                read.set_tag("SE",'NoFeatures')
            outsam.write(read)
        else:
            geneid = [x.value["geneid"] for x in overlapgene]
            read.set_tag("SG", "Ambiguity:[{0}]".format(",".join(geneid)))
            exonid = []
            for ii in overlapgene:
                overlapexon = getOverlapExon(read,flag,ii.value['exon'][chrom])
                if len(overlapexon) > 0:
                    exonid.append(ii.value["geneid"])
            if len(exonid) == 0:
                read.set_tag("SE",'NoFeatures')
            elif len(exonid) == 1:
                read.set_tag("SE",exonid[0])
            else:
                read.set_tag("SE","Ambiguity:[{0}]".format(",".join(exonid)))
            outsam.write(read)

if __name__ == "__main__":
    # sam/bam file, gtf file, outfile, flag
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])