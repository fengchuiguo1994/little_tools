import sys
import argparse
import itertools
import warnings
import traceback
import os.path

import pysam
import HTSeq

class UnknownChrom(Exception):
    pass

def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

def main(insamfile,gff_filename,outsamfile):
    com = ('M', '=', 'X')
    overlap_mode = "intersection-nonempty"
    multimapped_mode = "none"
    id_attribute = "gene_id"

    features_gene_sense = HTSeq.GenomicArrayOfSets("auto", True)
    features_exon_sense = HTSeq.GenomicArrayOfSets("auto", True)
    features_gene = HTSeq.GenomicArrayOfSets("auto", False)
    features_exon = HTSeq.GenomicArrayOfSets("auto", False)
    features_gene_antisense = HTSeq.GenomicArrayOfSets("auto", True)
    features_exon_antisense = HTSeq.GenomicArrayOfSets("auto", True)

    counts_gene_sense = {}
    counts_gene = {}
    counts_gene_antisense = {}
    counts_exon_sense = {}
    counts_exon = {}
    counts_exon_antisense = {}
    indexNumber = 0
    gff = HTSeq.GFF_Reader(gff_filename)
    for f in gff:
        if f.type == "gene":
            feature_id = f.attr[id_attribute]
            features_gene_sense[f.iv] += feature_id
            counts_gene_sense[f.attr[id_attribute]] = 0
            features_gene[f.iv] += feature_id
            counts_gene[f.attr[id_attribute]] = 0
            features_gene_antisense[f.iv] += feature_id
            counts_gene_antisense[f.attr[id_attribute]] = 0

        if f.type == "exon":
            feature_id = f.attr[id_attribute]
            features_exon_sense[f.iv] += feature_id
            counts_exon_sense[f.attr[id_attribute]] = 0
            features_exon[f.iv] += feature_id
            counts_exon[f.attr[id_attribute]] = 0
            features_exon_antisense[f.iv] += feature_id
            counts_exon_antisense[f.attr[id_attribute]] = 0

        indexNumber += 1
        if indexNumber % 100000 == 0:
            sys.stderr.write("%d GFF lines processed.\n" % indexNumber)
            sys.stderr.flush()
    if len(counts_gene) == 0:
        sys.exit("Warning: No features of type found.\n")

    if insamfile.endswith(".sam"):
        mysam = pysam.AlignmentFile(insamfile,'r')
        header = mysam.header
        mysam.close()
    else:
        mysam = pysam.AlignmentFile(insamfile,'rb')
        header = mysam.header
        mysam.close()

    if outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
        outsam.close()
    # else:
    #     outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    samoutfile = open(outsamfile, 'w+')

    if insamfile.endswith(".sam"):
        read_seq_file = HTSeq.SAM_Reader(insamfile)
    else:
        read_seq_file = HTSeq.BAM_Reader(insamfile)
    read_seq = iter(read_seq_file)

    indexNumber = 0
    for r in read_seq:
        if indexNumber > 0 and indexNumber % 100000 == 0:
            sys.stderr.write("{0} SAM alignment record processed.\n".format(indexNumber))
            sys.stderr.flush()
        indexNumber += 1

        if not r.aligned:
            r.optional_fields.append(('SG', "__not_aligned"))
            r.optional_fields.append(('SE', "__not_aligned"))
            r.optional_fields.append(('UG', "__not_aligned"))
            r.optional_fields.append(('UE', "__not_aligned"))
            r.optional_fields.append(('AG', "__not_aligned"))
            r.optional_fields.append(('AE', "__not_aligned"))
            samoutfile.write(r.get_sam_line() + "\n")
            continue
        # iv_seq_yes = (co.ref_iv for co in r.cigar if co.type in com and co.size > 0)
        # iv_seq_reverse = (invert_strand(co.ref_iv) for co in r.cigar if (co.type in com and co.size > 0))
        
        try:
            iv_seq_yes = [co.ref_iv for co in r.cigar if co.type in com and co.size > 0]
            iv_seq_reverse = [invert_strand(co.ref_iv) for co in r.cigar if (co.type in com and co.size > 0)]
            ##  sense gene ##
            fs = None
            for iv in iv_seq_yes:
                if iv.chrom not in features_gene_sense.chrom_vectors:
                    raise UnknownChrom
                for iv2, fs2 in features_gene_sense[iv].steps():
                    if (len(fs2) > 0):
                        if fs is None:
                            fs = fs2.copy()
                        else:
                            fs = fs.intersection(fs2)
            if fs is None or len(fs) == 0:
                r.optional_fields.append(('SG', "__no_feature"))
            elif len(fs) > 1:
                r.optional_fields.append(('SG', "__ambiguous[" + '+'.join(fs) + "]"))
            else:
                r.optional_fields.append(('SG', list(fs)[0]))

            ## sense exon ##
            fs = None
            for iv in iv_seq_yes:
                if iv.chrom not in features_exon_sense.chrom_vectors:
                    raise UnknownChrom
                for iv2, fs2 in features_exon_sense[iv].steps():
                    if (len(fs2) > 0):
                        if fs is None:
                            fs = fs2.copy()
                        else:
                            fs = fs.intersection(fs2)
            if fs is None or len(fs) == 0:
                r.optional_fields.append(('SE', "__no_feature"))
            elif len(fs) > 1:
                r.optional_fields.append(('SE', "__ambiguous[" + '+'.join(fs) + "]"))
            else:
                r.optional_fields.append(('SE', list(fs)[0]))

            ## unsense gene ##
            fs = None
            for iv in iv_seq_yes:
                if iv.chrom not in features_gene.chrom_vectors:
                    raise UnknownChrom
                for iv2, fs2 in features_gene[iv].steps():
                    if (len(fs2) > 0):
                        if fs is None:
                            fs = fs2.copy()
                        else:
                            fs = fs.intersection(fs2)
            if fs is None or len(fs) == 0:
                r.optional_fields.append(('UG', "__no_feature"))
            elif len(fs) > 1:
                r.optional_fields.append(('UG', "__ambiguous[" + '+'.join(fs) + "]"))
            else:
                r.optional_fields.append(('UG', list(fs)[0]))
            
            ## unsense exon ##
            fs = None
            for iv in iv_seq_yes:
                if iv.chrom not in features_exon.chrom_vectors:
                    raise UnknownChrom
                for iv2, fs2 in features_exon[iv].steps():
                    if (len(fs2) > 0):
                        if fs is None:
                            fs = fs2.copy()
                        else:
                            fs = fs.intersection(fs2)
            if fs is None or len(fs) == 0:
                r.optional_fields.append(('UE', "__no_feature"))
            elif len(fs) > 1:
                r.optional_fields.append(('UE', "__ambiguous[" + '+'.join(fs) + "]"))
            else:
                r.optional_fields.append(('UE', list(fs)[0]))
            
            ## antisense gene ##
            fs = None
            for iv in iv_seq_reverse:
                if iv.chrom not in features_gene_antisense.chrom_vectors:
                    raise UnknownChrom
                for iv2, fs2 in features_gene_antisense[iv].steps():
                    if (len(fs2) > 0):
                        if fs is None:
                            fs = fs2.copy()
                        else:
                            fs = fs.intersection(fs2)
            if fs is None or len(fs) == 0:
                r.optional_fields.append(('AG', "__no_feature"))
            elif len(fs) > 1:
                r.optional_fields.append(('AG', "__ambiguous[" + '+'.join(fs) + "]"))
            else:
                r.optional_fields.append(('AG', list(fs)[0]))

            ## antisense exon ##
            fs = None
            for iv in iv_seq_reverse:
                if iv.chrom not in features_exon_antisense.chrom_vectors:
                    raise UnknownChrom
                for iv2, fs2 in features_exon_antisense[iv].steps():
                    if (len(fs2) > 0):
                        if fs is None:
                            fs = fs2.copy()
                        else:
                            fs = fs.intersection(fs2)
            if fs is None or len(fs) == 0:
                r.optional_fields.append(('AE', "__no_feature"))
            elif len(fs) > 1:
                r.optional_fields.append(('AE', "__ambiguous[" + '+'.join(fs) + "]"))
            else:
                r.optional_fields.append(('AE', list(fs)[0]))
            
            samoutfile.write(r.get_sam_line() + "\n")

        except UnknownChrom:
            r.optional_fields.append(('SG', "__no_feature"))
            r.optional_fields.append(('SE', "__no_feature"))
            r.optional_fields.append(('UG', "__no_feature"))
            r.optional_fields.append(('UE', "__no_feature"))
            r.optional_fields.append(('AG', "__no_feature"))
            r.optional_fields.append(('AE', "__no_feature"))
            samoutfile.write(r.get_sam_line() + "\n")

    sys.stderr.write("{0} SAM reads processed.\n".format(indexNumber))
    sys.stderr.flush()
    samoutfile.close()


if __name__ == "__main__":
    # input file, gtf or gff file, out sam
    main(sys.argv[1],sys.argv[2],sys.argv[3])