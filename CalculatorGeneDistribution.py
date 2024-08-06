"""
针对bedpe文件，统计在全基因组上的分布
chr1 start1 end1 chr2 start2 end2 id score strand1 strand2 match1 match2 anchor gene
"""
import sys

def main(infile,outfile):
    # load file and store the dict
    genedict = {}
    genelocal = {}
    with open(infile,'r') as fin:
        for line in fin:
            tmp = line.strip().split()
            for i in tmp[13].split(";"):
                gene = i.split("-")[0]
                if gene not in genedict:
                    genedict[gene] = {}
                if tmp[0] not in genedict[gene]:
                    genedict[gene][tmp[0]] = 0
                genedict[gene][tmp[0]] += 1
                genelocal[gene] = tmp[3]

    # calculator and print output
    with open(outfile,'w') as fout:
        for gene in genedict.keys():
            tmpdict = genedict[gene]
            maxc = sumc = 0
            flag = None
            for chrom in tmpdict.keys():
                if tmpdict[chrom] > maxc:
                    maxc = tmpdict[chrom]
                    flag = chrom
                sumc += tmpdict[chrom]
            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene,genelocal[gene],chrom,maxc,sumc))

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2])