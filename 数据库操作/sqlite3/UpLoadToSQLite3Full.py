#!/usr/bin/python

import argparse
import sys
import os
import gzip
import re

import sqlite3
import shutil

precision = 4
header = "header"
databasedir = "database"
datadir = "data"
sizedir = "size"
chrbanddir = "chrband"
tablename = "table"
genomedir = "genome"

chrbandcolors = {'gpos1': '#FDFDFD', 'gpos2': '#FBFBFB', 'gpos3': '#F8F8F8', 'gpos4': '#F6F6F6', 'gpos5': '#F3F3F3', 'gpos6': '#F1F1F1', 'gpos7': '#EEEEEE', 'gpos8': '#ECECEC', 'gpos9': '#E9E9E9', 'gpos10': '#E6E6E6', 'gpos11': '#E4E4E4', 'gpos12': '#E1E1E1', 'gpos13': '#DFDFDF', 'gpos14': '#DCDCDC', 'gpos15': '#DADADA', 'gpos16': '#D7D7D7', 'gpos17': '#D4D4D4', 'gpos18': '#D2D2D2', 'gpos19': '#CFCFCF', 'gpos20': '#CDCDCD', 'gpos21': '#CACACA', 'gpos22': '#C8C8C8', 'gpos23': '#C5C5C5', 'gpos24': '#C3C3C3', 'gpos25': '#C0C0C0', 'gpos26': '#BDBDBD', 'gpos27': '#BBBBBB', 'gpos28': '#B8B8B8', 'gpos29': '#B6B6B6', 'gpos30': '#B3B3B3', 'gpos31': '#B1B1B1', 'gpos32': '#AEAEAE', 'gpos33': '#ACACAC', 'gpos34': '#A9A9A9', 'gpos35': '#A6A6A6', 'gpos36': '#A4A4A4', 'gpos37': '#A1A1A1', 'gpos38': '#9F9F9F', 'gpos39': '#9C9C9C', 'gpos40': '#9A9A9A', 'gpos41': '#979797', 'gpos42': '#949494', 'gpos43': '#929292', 'gpos44': '#8F8F8F', 'gpos45': '#8D8D8D', 'gpos46': '#8A8A8A', 'gpos47': '#888888', 'gpos48': '#858585', 'gpos49': '#838383', 'gpos50': '#808080', 'gpos51': '#7D7D7D', 'gpos52': '#7B7B7B', 'gpos53': '#787878', 'gpos54': '#767676', 'gpos55': '#737373', 'gpos56': '#717171', 'gpos57': '#6E6E6E', 'gpos58': '#6C6C6C', 'gpos59': '#696969', 'gpos60': '#666666', 'gpos61': '#646464', 'gpos62': '#616161', 'gpos63': '#5F5F5F', 'gpos64': '#5C5C5C', 'gpos65': '#5A5A5A', 'gpos66': '#575757', 'gpos67': '#545454', 'gpos68': '#525252', 'gpos69': '#4F4F4F', 'gpos70': '#4D4D4D', 'gpos71': '#4A4A4A', 'gpos72': '#484848', 'gpos73': '#454545', 'gpos74': '#434343', 'gpos75': '#404040', 'gpos76': '#3D3D3D', 'gpos77': '#3B3B3B', 'gpos78': '#383838', 'gpos79': '#363636', 'gpos80': '#333333', 'gpos81': '#313131', 'gpos82': '#2E2E2E', 'gpos83': '#2C2C2C', 'gpos84': '#292929', 'gpos85': '#262626', 'gpos86': '#242424', 'gpos87': '#212121', 'gpos88': '#1F1F1F', 'gpos89': '#1C1C1C', 'gpos90': '#1A1A1A', 'gpos91': '#171717', 'gpos92': '#141414', 'gpos93': '#121212', 'gpos94': '#0F0F0F', 'gpos95': '#0D0D0D', 'gpos96': '#0A0A0A', 'gpos97': '#080808', 'gpos98': '#050505', 'gpos99': '#030303', 'gpos100': '#000000', 'gneg': '#FFFFFF', 'acen': '#660033', 'gvar': '#660099', 'stalk': '#6600CC'}

# format:
# bigwig: binary
# hic: binary
# bedgraph: chrom start end value1
# ssbedgraph: chrom start end value1 value2 value3 ...
# ssbedgraph: chrom start end value1 value2
# bed (bed12): chrom start end name score strand thickstart thickend blockcount blocksizes blockstarts
# vcf: chrom start id ref alt qual filter(PASS/UNPASS) info
# indvcf: chrom start id ref alt qual filter(PASS/UNPASS) info
# bedpe: chrom1 start1 end1 chrom2 start2 end2 score/PETcount qvalue/pvalue
# scRNA: chrom start end barcode count
# scATAC: chrom start end barcode count
# scPET: chrom1 start1 end1 chrom2 start2 end2 barcode count
# mbed: chrom start end barcode count
# anno: chrom start end transcriptid score strand geneid type
# chromhmm: chrom start end chromstat score strand color
# size: # chrband: 

# meta/header: add group, assemble record

def traverse_directory(directory):
    refiles = []
    for root, dirs, files in os.walk(directory):
        for file_name in files:
            refiles.append(os.path.join(root, file_name))
    return(refiles)

def checkbed(mylist): # start1==start2  end1?end2
    return int(mylist[1]) >= int(mylist[2])

def compareChar(a1, a2):
    if a1 > a2:
        return 1
    elif a1 == a2:
        return 0
    else:
        return -1
    
def checkbedpe(mylist):
    result = compareChar(mylist[0], mylist[3])
    if result == 1:
        return True
    elif result == -1:
        return False
    else:
        return int(mylist[1]) >= int(mylist[2]) # start1==start2  end1?end2

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

def is_number(s):
    pattern = r'^-?\d+(\.\d+)?$'
    return bool(re.match(pattern, s))

import hicstraw
import pyBigWig
def detF(ifile):
    if ifile.endswith(".hic"):
        try:
            sys.stderr.write("check {0} file. ".format(ifile))
            sys.stderr.flush()
            hic_obj = hicstraw.HiCFile(ifile)
            return("hic")
        except:
            pass
        sys.exit("ERROR: {0} file is not a hic file".format(ifile))
    elif ifile.endswith((".bw", ".bigWig", ".BigWig", ".bigwig")):
        try:
            bb = pyBigWig.open(ifile)
            if bb.isBigWig():
                return("bw")
        except:
            pass
        sys.exit("ERROR: {0} file is not a bigwig file".format(ifile))
    elif ifile.endswith((".bedgraph", ".bedgraph.gz", ".bedGraph", ".bedGraph.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and is_number(tmp[3]):
            return("bedGraph")
        sys.exit("ERROR: {0} file is not a bedGraph file".format(ifile))
    elif ifile.endswith((".ssbedgraph", ".ssbedgraph.gz", ".ssbedGraph", ".ssbedGraph.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and is_number(tmp[3]) and is_number(tmp[4]):
            return("ssbedGraph")
        sys.exit("ERROR: {0} file is not a ssbedGraph file".format(ifile))
    elif ifile.endswith((".chromhmm", ".chromhmm.gz", ".chromHMM", ".chromHMM.gz")):
        return("chromhmm")
    elif ifile.endswith((".anno", ".anno.gz")):
        return("anno")
    elif ifile.endswith((".mbed", ".mbed.gz", ".Mbed", ".Mbed.gz")):
        return("mbed")
    elif ifile.endswith((".scPET", ".scPET.gz")):
        return("scPET")
    elif ifile.endswith((".scATAC", ".scATAC.gz")):
        return("scATAC")
    elif ifile.endswith((".scRNA", ".scRNA.gz")):
        return("scRNA")
    elif ifile.endswith((".bedpe", ".bedpe.gz", ".loop", ".loop.gz", ".curve", ".curve.gz")):
        return("bedpe")
    elif ifile.endswith((".indvcf",".indvcf.gz", ".indelvcf",".indelvcf.gz")):
        return("vcf")
    elif ifile.endswith((".vcf",".vcf.gz")):
        return("vcf")
    elif ifile.endswith((".bed", ".bed.gz", ".bed12", ".bed12.gz", ".bed6", ".bed6.gz")):
        return("bed")
    
    
    
    sys.exit("Please give a recognizable suffix")

def check(args, commands):
    fout = writeFile("PendingFileUpload.list")
    fout.write("#Source_path\tSource_type\tTarget_path\tTarget_format\tAssembly\n")
    if os.path.isfile(args.input): # file 
        basename =  os.path.basename(args.input)
        # file_name_without_extension, _ = os.path.splitext(basename.rstrip('.gz').rstrip('.GZ').rstrip('.gzip').rstrip('.GZIP'))
        file_name_without_extension = basename.rstrip('.gz').rstrip('.GZ').rstrip('.gzip').rstrip('.GZIP').rstrip('.Gzip')
        if args.format != "AUTO":
            if args.library == "":
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(args.input, args.format, file_name_without_extension, args.format, args.assembly))
            else:
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(args.input, args.format, os.path.join(args.library, file_name_without_extension), args.format, args.assembly))
        else:
            formatt = detF(args.input)
            if args.library == "":
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(args.input, formatt, file_name_without_extension, formatt, args.assembly))
            else:
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(args.input, formatt, os.path.join(args.library, file_name_without_extension), formatt, args.assembly))
    elif os.path.isdir(args.input): # directory
        for ii in traverse_directory(args.input):
            # file_name_without_extension, _ = os.path.splitext(ii)
            # file_name_without_extension, _ = os.path.splitext(ii.rstrip('.gz').rstrip('.GZ').rstrip('.gzip').rstrip('.GZIP'))
            file_name_without_extension = ii.rstrip('.gz').rstrip('.GZ').rstrip('.gzip').rstrip('.GZIP').rstrip('.Gzip')
            formatt = detF(ii)
            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(ii, formatt, file_name_without_extension, formatt, args.assembly))
    else:
        sys.exit("Not a valid path\n")
    fout.close()

def uploadbinary(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, "{0}".format(foldernumber))
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    shutil.copy(mylist[0], os.path.join(tmppath, "{0}.bin".format(foldernumber)))
    return 0

bedgraphformat = '''CREATE TABLE IF NOT EXISTS [{0}]
    (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
    chromosome     CHAR(100)         NOT NULL,
    start          INT               NOT NULL,
    end            INT               NOT NULL,
    value1         FLOAT             NOT NULL,
    value2         FLOAT             NULL);'''.format(tablename)
bedgraphformatindex = '''CREATE INDEX idx_{0} ON {0} (chromosome, start, end);'''.format(tablename)
def uploadbedgraph(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, "{0}".format(foldernumber))
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    tmppath = os.path.join(tmppath, "{0}.sqlite3".format(foldernumber))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    ###############################################################
    #           bedgraph (bulk RNA-Seq/ATAC-Seq...)               #
    ###############################################################
    c.execute(bedgraphformat)
    c.execute(bedgraphformatindex)
    conn.commit()
    fin = readFile(mylist[0])
    index = 0
    recordlist = []
    for line in fin:
        index += 1
        record = line.strip().split()
        if checkbed(record):
            sys.stderr.write("{0} line record is error!, please check\n".format(index))
            sys.stderr.flush()
            sys.exit("{0}\n".format("\t".join(record)))
        record[3] =  round(float(record[3]), precision)
        if mylist[3] == "ssbedGraph":
            record[4] =  round(float(record[4]), precision)
        recordlist.append(tuple(record))
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
            # print(recordlist[0])
            # print(mylist)
            if mylist[3] == "bedGraph":
                c.executemany("INSERT INTO [{0}] (chromosome, start, end, value1) VALUES (?,?,?,?)".format(tablename), recordlist)
            else:
                c.executemany("INSERT INTO [{0}] (chromosome, start, end, value1, value2) VALUES (?,?,?,?,?)".format(tablename), recordlist)
            conn.commit()
            recordlist = []
    if len(recordlist) > 0:
        # print(recordlist[0])
        # print(mylist)
        if mylist[3] == "bedGraph":
            c.executemany("INSERT INTO [{0}] (chromosome, start, end, value1) VALUES (?,?,?,?)".format(tablename), recordlist)
        else:
            c.executemany("INSERT INTO [{0}] (chromosome, start, end, value1, value2) VALUES (?,?,?,?,?)".format(tablename), recordlist)
        conn.commit()
    c.close()
    conn.close()
    return index

bed6format = '''CREATE TABLE IF NOT EXISTS [{0}]
    (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
    chromosome     CHAR(100)         NOT NULL,
    start          INT               NOT NULL,
    end            INT               NOT NULL,
    name           CHAR(100)         NOT NULL,
    score          INT               NOT NULL,
    strand         CHAR(1)           NOT NULL,
    color          CHAR(20));'''.format(tablename)
bed6formatindex = '''CREATE INDEX idx_{0} ON {0} (chromosome, start, end);'''.format(tablename)
def uploadbed6(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, foldernumber)
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    tmppath = os.path.join(tmppath, "{0}.sqlite3".format(foldernumber))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    ###############################################################
    #        bed6 (peak) /chromhmm (chrom stats) / chrband        #
    #                    scRNA / scATAC / mbed                    #
    ###############################################################
    c.execute(bed6format)
    c.execute(bed6formatindex)
    conn.commit()
    fin = readFile(mylist[0])
    index = 0
    recordlist = []
    for line in fin:
        index += 1
        record = line.strip().split()
        if checkbed(record):
            sys.stderr.write("{0} line record is error!, please check\n".format(index))
            sys.stderr.flush()
            sys.exit("{0}\n".format("\t".join(record)))
        if record[5] == "-" or record[5] == "0":
            record[5] == "-"
        elif record[5] == "+" or record[5] == "." or record[5] == "1":
            record[5] == "+"
        recordlist.append(tuple(record))
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
            if len(recordlist[0]) == 7:
                c.executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand, color) VALUES (?,?,?,?,?,?,?)".format(tablename), recordlist)
            elif len(recordlist[0]) == 6:
                c.executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand) VALUES (?,?,?,?,?,?)".format(tablename), recordlist)
            conn.commit()
            recordlist = []
    if len(recordlist) > 0:
        if len(recordlist[0]) == 7:
            c.executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand, color) VALUES (?,?,?,?,?,?,?)".format(tablename), recordlist)
        elif len(recordlist[0]) == 6:
            c.executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand) VALUES (?,?,?,?,?,?)".format(tablename), recordlist)
        conn.commit()
    c.close()
    conn.close()

vcfformat = '''CREATE TABLE IF NOT EXISTS [{0}]
    (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
    chromosome     CHAR(100)         NOT NULL,
    start          INT               NOT NULL,
    Mid            CHAR(100)         NOT NULL,
    ref            CHAR(100)         NOT NULL,
    alt            CHAR(100)         NOT NULL,
    qual           INT               NOT NULL);'''.format(tablename)
vcfformatindex = '''CREATE INDEX idx_{0} ON {0} (chromosome, start);'''.format(tablename)
def uploadvcf(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, foldernumber)
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    tmppath = os.path.join(tmppath, "{0}.sqlite3".format(foldernumber))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    ###############################################################
    #                           snp/indel                         #
    ###############################################################
    c.execute(vcfformat)
    c.execute(vcfformatindex)
    conn.commit()
    fin = readFile(mylist[0])
    index = 0
    recordlist = []
    for line in fin:
        index += 1
        record = line.strip().split()
        if checkbed(record):
            sys.stderr.write("{0} line record is error!, please check\n".format(index))
            sys.stderr.flush()
            sys.exit("{0}\n".format("\t".join(record)))
        recordlist.append(tuple(record))
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
            c.executemany("INSERT INTO [{0}] (chromosome, start, Mid, ref, alt, qual) VALUES (?,?,?,?,?,?)".format(tablename), recordlist)
            conn.commit()
            recordlist = []
    if len(recordlist) > 0:
        c.executemany("INSERT INTO [{0}] (chromosome, start, Mid, ref, alt, qual) VALUES (?,?,?,?,?,?)".format(tablename), recordlist)
        conn.commit()
    c.close()
    conn.close()

bedpeformat = '''CREATE TABLE IF NOT EXISTS [{0}]
    (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
    chromosome1     CHAR(100)         NOT NULL,
    start1          INT               NOT NULL,
    end1            INT               NOT NULL,
    chromosome2     CHAR(100)         NOT NULL,
    start2          INT               NOT NULL,
    end2            INT               NOT NULL,
    value           FLOAT             NOT NULL,
    barcode         CHAR(100),
    color           CHAR(10));'''.format(tablename)
bedpeformatindex = '''CREATE INDEX idx_{0} ON {0} (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(tablename)
def uploadbedpe(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, foldernumber)
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    tmppath = os.path.join(tmppath, "{0}.sqlite3".format(foldernumber))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    ###############################################################
    #                         loop/pet/scpet                      #
    ###############################################################
    c.execute(bedpeformat)
    c.execute(bedpeformatindex)
    conn.commit()
    fin = readFile(mylist[0])
    index = 0
    recordlist = []
    for line in fin:
        index += 1
        record = line.strip().split()
        if checkbedpe(record):
            sys.stderr.write("{0} line record is error!, please check\n".format(index))
            sys.stderr.flush()
            sys.exit("{0}\n".format("\t".join(record)))
        record[6] =  round(float(record[6]), precision)
        recordlist.append(tuple(record))
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
            if len(recordlist[0]) == 7:
                c.executemany("INSERT INTO [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value) VALUES (?,?,?,?,?,?,?)".format(tablename), recordlist)
            elif len(recordlist[0]) == 9:
                c.executemany("INSERT INTO [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value, barcode, color) VALUES (?,?,?,?,?,?,?,?,?)".format(tablename), recordlist)
            conn.commit()
            recordlist = []
    if len(recordlist) > 0:
        if len(recordlist[0]) == 7:
            c.executemany("INSERT INTO [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value) VALUES (?,?,?,?,?,?,?)".format(tablename), recordlist)
        elif len(recordlist[0]) == 9:
            c.executemany("INSERT INTO [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value, barcode, color) VALUES (?,?,?,?,?,?,?,?,?)".format(tablename), recordlist)
        conn.commit()
    c.close()
    conn.close()

annoformat = '''CREATE TABLE IF NOT EXISTS [{0}]
    (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
    chromosome     CHAR(100)         NOT NULL,
    start          INT               NOT NULL,
    end            INT               NOT NULL,
    geneID         CHAR(100)         NOT NULL,
    geneName       CHAR(100)         NOT NULL,
    strand         CHAR(1)           NOT NULL,
    type           CHAR(100)         NOT NULL);'''.format(tablename)
annoformatindex = '''CREATE INDEX idx_{0} ON {0} (chromosome, start, end);'''.format(tablename)
def uploadanno(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, foldernumber)
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    tmppath = os.path.join(tmppath, "{0}.sqlite3".format(foldernumber))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    ###############################################################
    #                            anno                             #
    ###############################################################
    c.execute(annoformat)
    c.execute(annoformatindex)
    conn.commit()
    fin = readFile(mylist[0])
    index = 0
    recordlist = []
    for line in fin:
        index += 1
        record = line.strip().split()
        if checkbed(record):
            sys.stderr.write("{0} line record is error!, please check\n".format(index))
            sys.stderr.flush()
            sys.exit("{0}\n".format("\t".join(record)))
        if record[5] == "-" or record[5] == "0":
            record[5] == "-"
        elif record[5] == "+" or record[5] == "." or record[5] == "1":
            record[5] == "+"
        recordlist.append(tuple(record))
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
            c.executemany("INSERT INTO [{0}] (chromosome, start, end, geneID, geneName, strand, type) VALUES (?,?,?,?,?,?,?)".format(tablename), recordlist)
            conn.commit()
            recordlist = []
    if len(recordlist) > 0:
        c.executemany("INSERT INTO [{0}] (chromosome, start, end, geneID, geneName, strand, type) VALUES (?,?,?,?,?,?,?)".format(tablename), recordlist)
        conn.commit()
    c.close()
    conn.close()

def upload(args, commands):
    datapath = os.path.join(databasedir, datadir)
    if not os.path.exists(datapath):
        os.makedirs(datapath)
    all_items = os.listdir(datapath)
    folder_names = [int(item) for item in all_items if os.path.isdir(os.path.join(datapath, item))]
    if len(folder_names) == 0:
        curdat = 0
    else:
        curdat = max(folder_names)
    fin = readFile(args.input)
    fout = writeFile("{0}.upload".format(args.input))
    conn = sqlite3.connect("{0}.sqlite3".format(datapath))
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS [{0}]
        (ID           INTEGER PRIMARY KEY   AUTOINCREMENT,
        sample        CHAR(100)             NOT NULL,
        format        CHAR(20)              NOT NULL,
        fileName      CHAR(100)             NOT NULL,
        assembly      CHAR(100)             NOT NULL,
        nrecords      INTEGER               NOT NULL,
        dataID        INTEGER               NOT NULL,
        groupID       CHAR(100)             NULL);'''.format(tablename))
    conn.commit()
    for line in fin:
        if line.startswith("#"):
            fout.write("{0}\tid\n".format(line.strip()))
            continue
        curdat += 1
        tmp = line.strip().split()
        if tmp[3] == "bw" or tmp[3] == "hic":
            nrecords = uploadbinary(tmp, curdat, datapath)
        elif tmp[3] == "bedGraph" or tmp[3] == "ssbedGraph":
            nrecords = uploadbedgraph(tmp, curdat, datapath)
        elif tmp[3] == "bed6" or tmp[3] == "scRNA" or tmp[3] == "scATAC" or tmp[3] == "mbed" or tmp[3] == "chromhmm":
            uploadbed6(tmp, curdat, datapath)
        elif tmp[3] == "vcf":
            uploadvcf(tmp, curdat, datapath)
        elif tmp[3] == "bedpe" or tmp[3] == "scPET":
            uploadbedpe(tmp, curdat, datapath)
        elif tmp[3] == "anno":
            uploadanno(tmp, curdat, datapath)
        else:
            sys.exit("file format error")
        c.execute('INSERT INTO [{0}] (sample, format, fileName, assembly, nrecords, dataID) VALUES (?,?,?,?,?,?)'.format(tablename), (tmp[2], tmp[3], tmp[0], tmp[4], nrecords, curdat))
        conn.commit()
        fout.write("{0}\t{1}\n".format(line.strip(), curdat))
        fout.flush()
    fin.close()
    fout.close()
    c.close()
    conn.close()

def addGenome(args, commands):
    if not os.path.exists(databasedir):
        os.makedirs(databasedir)
    datapath = os.path.join(databasedir, datadir)
    if not os.path.exists(datapath):
        os.makedirs(datapath)
    genomepath = os.path.join(databasedir, genomedir)
    if not os.path.exists(genomepath):
        os.makedirs(genomepath)
    sizepath = os.path.join(genomepath, sizedir)
    if not os.path.exists(sizepath):
        os.makedirs(sizepath)
    chrbandpath = os.path.join(genomepath, chrbanddir)
    if not os.path.exists(chrbandpath):
        os.makedirs(chrbandpath)
    if os.path.exists(os.path.join(sizepath, "{0}.size".format(args.assembly))):
        sys.exit("ERROR: {0} is already exists, please change another names/assembly".format(args.assembly))
    
    fin = readFile(args.size)
    fout = writeFile(os.path.join(sizepath, "{0}.size".format(args.assembly)))
    for line in fin:
        if line.strip() == "":
            continue
        tmp = line.strip().split()
        fout.write("{0}\t{1}\n".format(tmp[0], tmp[1]))
    fin.close()
    fout.close()

    if args.chrband != "":
        fin = readFile(args.chrband)
        fout = writeFile(os.path.join(chrbandpath, "{0}.chrband".format(args.assembly)))
        for line in fin:
            if line.strip() == "":
                continue
            tmp = line.strip().split()
            fout.write("{0}\t{1}\n".format("\t".join(tmp[0:5]), chrbandcolors[tmp[4]]))
        fin.close()
        fout.close()
    else:
        fin = readFile(args.size)
        fout = writeFile(os.path.join(chrbandpath, "{0}.chrband".format(args.assembly)))
        for line in fin:
            if line.strip() == "":
                continue
            tmp = line.strip().split()
            fout.write("{0}\t0\t{1}\tp\tgneg\t{2}\n".format(tmp[0],tmp[1], chrbandcolors["gneg"]))
        fin.close()
        fout.close()

    conn = sqlite3.connect("{0}.sqlite3".format(os.path.join(databasedir, genomedir)))
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS [{0}]
        (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
        assembly       CHAR(100)             NOT NULL);'''.format(tablename))
    conn.commit()
    c.execute('INSERT INTO [{0}] (assembly) VALUES (?)'.format(tablename), (args.assembly,))
    conn.commit()
    c.close()
    conn.close()

def rmGenome(args, commands):
    genomepath = os.path.join(databasedir, genomedir)
    sizepath = os.path.join(genomepath, sizedir)
    if os.path.exists(os.path.join(sizepath, "{0}.size".format(args.assembly))):
        os.remove(os.path.join(sizepath, "{0}.size".format(args.assembly)))
    else:
        sys.exit("ERROR: {0} is not exists".format(args.assembly))
    chrbandpath = os.path.join(genomepath, chrbanddir)
    if os.path.exists(os.path.join(chrbandpath, "{0}.chrband".format(args.assembly))):
        os.remove(os.path.join(chrbandpath, "{0}.chrband".format(args.assembly)))
    else:
        sys.exit("ERROR: {0} is not exists".format(args.assembly))

    conn = sqlite3.connect("{0}.sqlite3".format(os.path.join(databasedir, genomedir)))
    c = conn.cursor()
    c.execute('DELETE FROM [{0}] WHERE assembly = ?'.format(tablename), (args.assembly,))
    conn.commit()
    c.close()
    conn.close()
    print("{0} remove successful".format(args.assembly))

def getargs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Program: UpLoadToSQLite3Full
Version: 1.0.0
Usage: UpLoadToSQLite3Full <command> [options]
Commands:
    filecheck   check the upload file, and generate the upload information.
    remove      remove table
==========================================================================

""")
    parser.add_argument('--version', action='version',version='%(prog)s {0} by xyhuang'.format('1.0.0'))
    subparsers = parser.add_subparsers(title='subcommands',description='valid subcommands',help='') 
    
    ## add assembly
    parser_add = subparsers.add_parser('addGenome',help='add a new assembly')
    parser_add.add_argument("-a", "--assembly", required=True, help = "assembly name, for example: mm10, hg38, hg19")
    parser_add.add_argument("-s", "--size", required=True, help="the genome size file")
    parser_add.add_argument("-c", "--chrband", required=False, default="" , help="the genome chrband file")
    parser_add.set_defaults(func=addGenome)

    ## remove assembly
    parser_add = subparsers.add_parser('rmGenome',help='remove assembly')
    parser_add.add_argument("-a", "--assembly", required=True, help = "assembly name, for example: mm10, hg38, hg19")
    parser_add.set_defaults(func=rmGenome)

    ## check upload file format
    parser_add = subparsers.add_parser('filecheck',help='check the upload file, and generate the upload information')
    parser_add.add_argument("-i", "--input", required=True, help="the input file/directory")
    parser_add.add_argument("-f", "--format", required=False, default="AUTO", choices=["AUTO", "bigwig", "hic", "bedgraph", "ssbedgraph", "bed", "vcf", "bedpe", "scRNA", "scATAC", "scPET", "mbed", "anno", "chromhmm"] ,help="Format of input file. If set to AUTO, the program will determine the file format based on the filename suffix. If the inputdirectory parameter is set, it defaults to AUTO [default:AUTO]")
    parser_add.add_argument("-a", "--assembly", required=True, help = "assembly, for example: mm10, hg38, hg19")
    parser_add.add_argument("-l", "--library", required=False, default="", help = "library/directory, the library in web. It works when the input is a file, and is ignored if it is a directory.")
    parser_add.set_defaults(func=check)

    ## upload file
    parser_add = subparsers.add_parser('upload',help='upload file')
    parser_add.add_argument("-i", "--input", required=True, help="the upload information (generate by filecheck)")
    parser_add.set_defaults(func=upload)

    commands = sys.argv[1:]
    if ((not commands) or ((commands[0] in ['add', 'remove']) and len(commands) == 1)):
        commands.append('-h')
    args = parser.parse_args(commands)
    return args, commands

if __name__ == "__main__":
    args, commands = getargs()
    args.func(args, commands)
