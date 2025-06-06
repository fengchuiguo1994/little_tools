import argparse
import sys
import pysam

def readSam(insamfile):
    if insamfile.endswith((".bam",".sam.gz")):
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

def detPM(mystr):
    if mystr == "noallele":
        return "noallele"
    mylist = mystr.split(";")
    if mylist.count("P") >= mylist.count("M") * 4: # 80%
        return "P"
    elif mylist.count("M") >= mylist.count("P") * 4:
        return "M"
    else:
        return "amb"

def build(read):
    a = pysam.AlignedSegment()
    a.query_name = read.query_name
    a.query_sequence = read.query_sequence
    a.flag = read.flag
    a.reference_start = read.reference_start
    a.mapping_quality = read.mapping_quality
    a.cigar = read.cigar
    a.next_reference_start = read.next_reference_start
    a.template_length = read.template_length
    a.query_qualities = read.query_qualities
    a.tags = read.tags
    return a

parser = argparse.ArgumentParser(description = "  ")
parser.add_argument("-i","--bamfile",required=True,help=" the input bam file ")
parser.add_argument("-o","--outputprefix",required=True,help=" output prefix ")
args = parser.parse_args()

if __name__ == "__main__":
    fin = readSam(args.bamfile)
    header = fin.header
    aa = header.as_dict()
    bb = aa["SQ"]
    cc = []
    for index,x in enumerate(bb):
        tmp = {}
        tmp["SN"] = "{0}-M".format(bb[index]["SN"])
        tmp["LN"] = bb[index]["LN"]
        cc.append(tmp)
        tmp = {}
        tmp["SN"] = "{0}-P".format(bb[index]["SN"])
        tmp["LN"] = bb[index]["LN"]
        cc.append(tmp)
        tmp = {}
        tmp["SN"] = "{0}".format(bb[index]["SN"])
        tmp["LN"] = bb[index]["LN"]
        cc.append(tmp)
    aa["SQ"] = cc

    # print(aa)
    filehash = {}
    filehash["flt"] = writeSam("{0}.flt.bam".format(args.outputprefix),aa)
    # print(filehash["flt"].header)
    stat = {}
    for read1 in fin:
        read2 = next(fin)
        flag1 = detPM(read1.get_tag("AL"))
        flag2 = detPM(read2.get_tag("AL"))
        if flag1 == "noallele" and flag2 == "noallele":
            a = build(read1)
            b = build(read2)
            a.reference_id = filehash["flt"].get_tid("{0}".format(read1.reference_name))
            a.next_reference_id = filehash["flt"].get_tid("{0}".format(read2.reference_name))
            # if read1.query_name == "457135029:24154:8519809:58786689:1433:29306:34147":
            #     print(read1.reference_name)
            #     print(filehash["flt"].get_reference_name(a.reference_id))
            b.reference_id = filehash["flt"].get_tid("{0}".format(read2.reference_name))
            b.next_reference_id = filehash["flt"].get_tid("{0}".format(read1.reference_name))
            filehash["flt"].write(a)
            filehash["flt"].write(b)
            if "noallele" not in stat:
                stat["noallele"] = 0
            stat["noallele"] += 1
        elif flag1 == "amb" or flag2 == "amb":
            a = build(read1)
            b = build(read2)
            a.reference_id = filehash["flt"].get_tid("{0}".format(read1.reference_name))
            a.next_reference_id = filehash["flt"].get_tid("{0}".format(read2.reference_name))
            b.reference_id = filehash["flt"].get_tid("{0}".format(read2.reference_name))
            b.next_reference_id = filehash["flt"].get_tid("{0}".format(read1.reference_name))
            filehash["flt"].write(a)
            filehash["flt"].write(b)
            if "amb" not in stat:
                stat["amb"] = 0
            stat["amb"] += 1
        else:
            a = build(read1)
            b = build(read2)
            if flag1 == "noallele":
                # print(a)
                # print(b)
                a.reference_id = filehash["flt"].get_tid("{0}".format(read1.reference_name))
                a.next_reference_id = filehash["flt"].get_tid("{0}-{1}".format(read2.reference_name,flag2))
                b.reference_id = filehash["flt"].get_tid("{0}-{1}".format(read2.reference_name,flag2))
                b.next_reference_id = filehash["flt"].get_tid("{0}".format(read1.reference_name))
                # print(a)
                # print(b)
            elif flag2 == "noallele":
                # print(a)
                # print(b)
                a.reference_id = filehash["flt"].get_tid("{0}-{1}".format(read1.reference_name,flag1))
                a.next_reference_id = filehash["flt"].get_tid("{0}".format(read2.reference_name))
                b.reference_id = filehash["flt"].get_tid("{0}".format(read2.reference_name))
                b.next_reference_id = filehash["flt"].get_tid("{0}-{1}".format(read1.reference_name,flag1))
                # print(a)
                # print(b)
            else:
                a.reference_id = filehash["flt"].get_tid("{0}-{1}".format(read1.reference_name,flag1))
                a.next_reference_id = filehash["flt"].get_tid("{0}-{1}".format(read2.reference_name,flag2))
                b.reference_id = filehash["flt"].get_tid("{0}-{1}".format(read2.reference_name,flag2))
                b.next_reference_id = filehash["flt"].get_tid("{0}-{1}".format(read1.reference_name,flag1))
            filehash["flt"].write(a)
            filehash["flt"].write(b)
    for k in sorted(filehash.keys()):
        filehash[k].close()
    fin.close()
    for i in stat:
        print("{0}\t{1}".format(i,stat[i]))