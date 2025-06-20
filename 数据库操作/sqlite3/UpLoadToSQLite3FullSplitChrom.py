#!/usr/bin/python

import argparse
import sys
import os
import gzip
import re
import random

import sqlite3
import shutil

precision = 4
header = "header"
databasedir = "database"
datadir = "data"
cifdatadir = "cifdata"
sizedir = "size"
chrbanddir = "chrband"
tablename = "table"
genomedir = "genome"
nrecord = 10000

chrbandcolors = {'gpos1': '#FDFDFD', 'gpos2': '#FBFBFB', 'gpos3': '#F8F8F8', 'gpos4': '#F6F6F6', 'gpos5': '#F3F3F3', 'gpos6': '#F1F1F1', 'gpos7': '#EEEEEE', 'gpos8': '#ECECEC', 'gpos9': '#E9E9E9', 'gpos10': '#E6E6E6', 'gpos11': '#E4E4E4', 'gpos12': '#E1E1E1', 'gpos13': '#DFDFDF', 'gpos14': '#DCDCDC', 'gpos15': '#DADADA', 'gpos16': '#D7D7D7', 'gpos17': '#D4D4D4', 'gpos18': '#D2D2D2', 'gpos19': '#CFCFCF', 'gpos20': '#CDCDCD', 'gpos21': '#CACACA', 'gpos22': '#C8C8C8', 'gpos23': '#C5C5C5', 'gpos24': '#C3C3C3', 'gpos25': '#C0C0C0', 'gpos26': '#BDBDBD', 'gpos27': '#BBBBBB', 'gpos28': '#B8B8B8', 'gpos29': '#B6B6B6', 'gpos30': '#B3B3B3', 'gpos31': '#B1B1B1', 'gpos32': '#AEAEAE', 'gpos33': '#ACACAC', 'gpos34': '#A9A9A9', 'gpos35': '#A6A6A6', 'gpos36': '#A4A4A4', 'gpos37': '#A1A1A1', 'gpos38': '#9F9F9F', 'gpos39': '#9C9C9C', 'gpos40': '#9A9A9A', 'gpos41': '#979797', 'gpos42': '#949494', 'gpos43': '#929292', 'gpos44': '#8F8F8F', 'gpos45': '#8D8D8D', 'gpos46': '#8A8A8A', 'gpos47': '#888888', 'gpos48': '#858585', 'gpos49': '#838383', 'gpos50': '#808080', 'gpos51': '#7D7D7D', 'gpos52': '#7B7B7B', 'gpos53': '#787878', 'gpos54': '#767676', 'gpos55': '#737373', 'gpos56': '#717171', 'gpos57': '#6E6E6E', 'gpos58': '#6C6C6C', 'gpos59': '#696969', 'gpos60': '#666666', 'gpos61': '#646464', 'gpos62': '#616161', 'gpos63': '#5F5F5F', 'gpos64': '#5C5C5C', 'gpos65': '#5A5A5A', 'gpos66': '#575757', 'gpos67': '#545454', 'gpos68': '#525252', 'gpos69': '#4F4F4F', 'gpos70': '#4D4D4D', 'gpos71': '#4A4A4A', 'gpos72': '#484848', 'gpos73': '#454545', 'gpos74': '#434343', 'gpos75': '#404040', 'gpos76': '#3D3D3D', 'gpos77': '#3B3B3B', 'gpos78': '#383838', 'gpos79': '#363636', 'gpos80': '#333333', 'gpos81': '#313131', 'gpos82': '#2E2E2E', 'gpos83': '#2C2C2C', 'gpos84': '#292929', 'gpos85': '#262626', 'gpos86': '#242424', 'gpos87': '#212121', 'gpos88': '#1F1F1F', 'gpos89': '#1C1C1C', 'gpos90': '#1A1A1A', 'gpos91': '#171717', 'gpos92': '#141414', 'gpos93': '#121212', 'gpos94': '#0F0F0F', 'gpos95': '#0D0D0D', 'gpos96': '#0A0A0A', 'gpos97': '#080808', 'gpos98': '#050505', 'gpos99': '#030303', 'gpos100': '#000000', 'gneg': '#FFFFFF', 'acen': '#660033', 'gvar': '#660099', 'stalk': '#6600CC'}

# format:
# bigwig: binary
# hic: binary
# mcool: binary
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

# import hicstraw
# import pyBigWig
# import cooler
# import h5py
def detF(ifile):
    if ifile.endswith(".hic"):
        # try:
        #     sys.stderr.write("check {0} file. ".format(ifile))
        #     sys.stderr.flush()
        #     hic_obj = hicstraw.HiCFile(ifile)
        #     return("hic")
        # except:
        #     pass
        # sys.exit("ERROR: {0} file is not a hic file".format(ifile))
        return("hic")
    elif ifile.endswith((".bw", ".bigWig", ".BigWig", ".bigwig")):
        # try:
        #     bb = pyBigWig.open(ifile)
        #     if bb.isBigWig():
        #         return("bw")
        # except:
        #     pass
        # sys.exit("ERROR: {0} file is not a bigwig file".format(ifile))
        return("bw")
    elif ifile.endswith(".mcool"):
        # try:
        #     bb = h5py.File(ifile, 'r')
        #     resolutions = list(bb['resolutions'].keys())
        #     Lib = cooler.Cooler("all.K562.hg38.mcool::resolutions/{0}".format(resolutions[0]))
        #     return("mcool")
        # except:
        #     pass
        # sys.exit("ERROR: {0} file is not a mcool file".format(ifile))
        return("mcool")
    elif ifile.endswith(".cif"):
        return("cif")
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
    elif ifile.endswith((".chromhmm", ".chromhmm.gz", ".chromHMM", ".chromHMM.gz", ".chromhmm.bed", ".chromhmm.bed.gz", ".chromHMM.bed", ".chromHMM.bed.gz", ".Chromhmm", ".Chromhmm.gz", ".ChromHMM", ".ChromHMM.gz", ".Chromhmm.bed", ".Chromhmm.bed.gz", ".ChromHMM.bed", ".ChromHMM.bed.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and len(tmp) >= 7:
            return("chromhmm")
        sys.exit("ERROR: {0} file is not a chromhmm file".format(ifile))
    elif ifile.endswith((".anno", ".anno.gz", ".anno.txt", ".anno.txt.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and len(tmp) >= 8:
            return("anno")
        sys.exit("ERROR: {0} file is not a anno file".format(ifile))
    elif ifile.endswith((".mbed", ".mbed.gz", ".Mbed", ".Mbed.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and len(tmp) >= 5:
            return("mbed")
        sys.exit("ERROR: {0} file is not a mbed file".format(ifile))
    elif ifile.endswith((".scATAC", ".scATAC.gz", ".scATAC.bed", ".scATAC.bed.gz", ".scATAC.txt", ".scATAC.txt.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and len(tmp) >= 5:
            return("scATAC")
        sys.exit("ERROR: {0} file is not a scATAC file".format(ifile))
    elif ifile.endswith((".scRNA", ".scRNA.gz", ".scRNA.bed", ".scRNA.bed.gz", ".scRNA.txt", ".scRNA.txt.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and len(tmp) >= 5:
            return("scRNA")
        sys.exit("ERROR: {0} file is not a scRNA file".format(ifile))
    elif ifile.endswith((".vcf",".vcf.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and len(tmp) >= 5:
            return("vcf")
        sys.exit("ERROR: {0} file is not a vcf file".format(ifile))
    elif ifile.endswith((".scPET", ".scPET.gz", ".scPET.bedpe", ".scPET.bedpe.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and tmp[4].isdigit() and tmp[5].isdigit() and len(tmp) >= 8:
            return("scPET")
        sys.exit("ERROR: {0} file is not a scPET file".format(ifile))
    elif ifile.endswith((".bedpe", ".bedpe.gz", ".loop", ".loop.gz", ".curve", ".curve.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and tmp[4].isdigit() and tmp[5].isdigit() and len(tmp) >= 7:
            return("bedpe")
        sys.exit("ERROR: {0} file is not a bedpe file".format(ifile))
    elif ifile.endswith((".bed", ".bed.gz", ".bed12", ".bed12.gz", ".bed6", ".bed6.gz")):
        fin = readFile(ifile)
        tmp = fin.readline().strip().split()
        if tmp[1].isdigit() and tmp[2].isdigit() and len(tmp) >= 6:
            return("bed6")
        sys.exit("ERROR: {0} file is not a bed6 file".format(ifile))
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
    if mylist[3] == "cif":
        tmppath = os.path.join(foldpath, "../{0}/{1}".format(cifdatadir, foldernumber))
    else:
        tmppath = os.path.join(foldpath, "{0}".format(foldernumber))
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    shutil.copy(mylist[0], os.path.join(tmppath, "{0}.bin".format(foldernumber)))
    return 0

def conndbbedgraph(tmppath, foldernumber, tablename, chromosome):
    tmppath = os.path.join(tmppath, "{0}.{1}.sqlite3".format(foldernumber, chromosome))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    bedgraphformat = '''CREATE TABLE IF NOT EXISTS [{0}]
        (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
        chromosome     CHAR(100)         NOT NULL,
        start          INT               NOT NULL,
        end            INT               NOT NULL,
        value1         FLOAT             NOT NULL,
        value2         FLOAT             NULL);'''.format(tablename)
    bedgraphformatindex = '''CREATE INDEX [idx_{0}] ON [{0}] (chromosome, start, end, value1, value2);'''.format(tablename)
    c.execute(bedgraphformat)
    c.execute(bedgraphformatindex)
    conn.commit()
    return conn, c

def uploadbedgraph(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, "{0}".format(foldernumber))
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    print(tmppath)
    ###############################################################
    #                     bedgraph / ssbedgraph                   #
    ###############################################################
    fin = readFile(mylist[0])
    recordhash = {}
    dbhash = {}
    index = 0
    for line in fin:
        record = line.strip().split()
        if checkbed(record):
            sys.stderr.write("{0} line record is error!, please check\n".format(index))
            sys.stderr.flush()
            sys.exit("{0}\n".format("\t".join(record)))
        record[3] =  round(float(record[3]), precision)
        if mylist[3] == "ssbedGraph":
            record[4] =  round(float(record[4]), precision)
        if record[0] not in recordhash:
            recordhash[record[0]] = []
            conn, c = conndbbedgraph(tmppath, foldernumber, tablename, record[0])
            dbhash[record[0]] = {}
            dbhash[record[0]]["conn"] = conn
            dbhash[record[0]]["c"] = c
        recordhash[record[0]].append(tuple(record))
        if len(recordhash[record[0]]) == nrecord:
            if mylist[3] == "bedGraph":
                dbhash[record[0]]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, value1) VALUES (?,?,?,?)".format(tablename), recordhash[record[0]])
            else:
                dbhash[record[0]]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, value1, value2) VALUES (?,?,?,?,?)".format(tablename), recordhash[record[0]])
            dbhash[record[0]]["conn"].commit()
            recordhash[record[0]] = []
        index += 1
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
    for k in recordhash.keys():
        if len(recordhash[k]) > 0:
            if mylist[3] == "bedGraph":
                dbhash[k]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, value1) VALUES (?,?,?,?)".format(tablename), recordhash[k])
            else:
                dbhash[k]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, value1, value2) VALUES (?,?,?,?,?)".format(tablename), recordhash[k])
            dbhash[k]["conn"].commit()
        dbhash[k]["c"].close()
        dbhash[k]["conn"].close()
    return index

def conndbbed(tmppath, foldernumber, tablename, chromosome):
    tmppath = os.path.join(tmppath, "{0}.{1}.sqlite3".format(foldernumber, chromosome))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    bed6format = '''CREATE TABLE IF NOT EXISTS [{0}]
        (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
        chromosome     CHAR(100)         NOT NULL,
        start          INT               NOT NULL,
        end            INT               NOT NULL,
        name           CHAR(100)         NOT NULL,
        score          INT               NOT NULL,
        strand         CHAR(1)           NOT NULL,
        value1         CHAR(20),
        value2         CHAR(20),
        Txstart        INT,
        Txend          INT);'''.format(tablename)
    bed6formatindex = '''CREATE INDEX [idx_{0}] ON [{0}] (chromosome, start, end, name, score, strand, value1, value2, Txstart, Txend);'''.format(tablename)
    c.execute(bed6format)
    c.execute(bed6formatindex)
    conn.commit()
    return conn, c

def uploadbed6(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, "{0}".format(foldernumber))
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    print(tmppath)
    ###############################################################
    #             bed6 / anno /chromhmm (chrom stats)             #
    #                    scRNA / scATAC / mbed                    #
    ###############################################################
    fin = readFile(mylist[0])
    recordhash = {}
    dbhash = {}
    index = 0
    barcodedict = {}
    for line in fin:
        record = line.strip().split()
        if checkbed(record):
            sys.stderr.write("{0} line record is error!, please check\n".format(index))
            sys.stderr.flush()
            sys.exit("{0}\n".format("\t".join(record)))
        if len(record) < 5:
            record.append(0)
        if len(record) < 6:
            record.append("+")
        if record[5] == "-" or record[5] == "0":
            record[5] == "-"
        if record[5] == "+" or record[5] == "." or record[5] == "1":
            record[5] == "+"
        if record[0] not in recordhash:
            recordhash[record[0]] = []
            conn, c = conndbbed(tmppath, foldernumber, tablename, record[0])
            dbhash[record[0]] = {}
            dbhash[record[0]]["conn"] = conn
            dbhash[record[0]]["c"] = c
        recordhash[record[0]].append(tuple(record))
        if mylist[3] == "mbed" or mylist[3] == "scATAC" or mylist[3] == "scRNA":
            if record[3] not in barcodedict:
                barcodedict[record[3]] = 0
            barcodedict[record[3]] += 1
        if len(recordhash[record[0]]) == nrecord:
            if mylist[3] == "chromhmm":
                dbhash[record[0]]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand, value1) VALUES (?,?,?,?,?,?,?)".format(tablename), recordhash[record[0]])
            elif mylist[3] == "anno":
                dbhash[record[0]]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand, value1, value2, Txstart, Txend) VALUES (?,?,?,?,?,?,?,?,?,?)".format(tablename), recordhash[record[0]])
            elif mylist[3] == "mbed" or mylist[3] == "scATAC" or mylist[3] == "scRNA" or mylist[3] == "bed6":
                dbhash[record[0]]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand) VALUES (?,?,?,?,?,?)".format(tablename), recordhash[record[0]])
            dbhash[record[0]]["conn"].commit()
            recordhash[record[0]] = []
        index += 1
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
    for k in recordhash.keys():
        if len(recordhash[k]) > 0:
            if mylist[3] == "chromhmm":
                dbhash[k]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand, value1) VALUES (?,?,?,?,?,?,?)".format(tablename), recordhash[k])
            elif mylist[3] == "anno":
                dbhash[k]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand, value1, value2, Txstart, Txend) VALUES (?,?,?,?,?,?,?,?,?,?)".format(tablename), recordhash[k])
            elif mylist[3] == "mbed" or mylist[3] == "scATAC" or mylist[3] == "scRNA" or mylist[3] == "bed6":
                dbhash[k]["c"].executemany("INSERT INTO [{0}] (chromosome, start, end, name, score, strand) VALUES (?,?,?,?,?,?)".format(tablename), recordhash[k])
            dbhash[k]["conn"].commit()
        dbhash[k]["c"].close()
        dbhash[k]["conn"].close()
    if mylist[3] == "mbed" or mylist[3] == "scATAC" or mylist[3] == "scRNA":
        connttt = sqlite3.connect(os.path.join(tmppath, "{0}.barcodeindex.sqlite3".format(foldernumber)))
        cttt = connttt.cursor()
        cttt.execute('''CREATE TABLE IF NOT EXISTS [{0}]
        (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
        barcode        CHAR(100)         NOT NULL,
        count          INT               NOT NULL);'''.format(tablename))
        connttt.commit()
        tmplist = []
        barcodelist = list(barcodedict.keys())
        # print(barcodelist[0:10])
        random.shuffle(barcodelist)
        # print(barcodelist[0:10])
        nindex = 0
        for k in barcodelist:
            tmplist.append((k, barcodedict[k]))
            nindex += 1
            if nindex % nrecord == 0:
                cttt.executemany('''INSERT INTO [{0}] (barcode, count) VALUES (?,?)'''.format(tablename), tmplist)
                tmplist = []
        if len(tmplist) > 0:
            cttt.executemany('''INSERT INTO [{0}] (barcode, count) VALUES (?,?)'''.format(tablename), tmplist)
        connttt.commit()
        cttt.close()
        connttt.close()
    sys.stdout.flush()
    return index

def conndbvcf(tmppath, foldernumber, tablename, chromosome):
    tmppath = os.path.join(tmppath, "{0}.{1}.sqlite3".format(foldernumber, chromosome))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    vcfformat = '''CREATE TABLE IF NOT EXISTS [{0}]
    (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
    chromosome     CHAR(100)         NOT NULL,
    start          INT               NOT NULL,
    Mid            CHAR(100)         NOT NULL,
    ref            CHAR(100)         NOT NULL,
    alt            CHAR(100)         NOT NULL);'''.format(tablename)
    vcfformatindex = '''CREATE INDEX [idx_{0}] ON [{0}] (chromosome, start, Mid, ref, alt);'''.format(tablename)
    c.execute(vcfformat)
    c.execute(vcfformatindex)
    conn.commit()
    return conn, c

def uploadvcf(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, "{0}".format(foldernumber))
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    print(tmppath)
    ###############################################################
    #                           snp/indel                         #
    ###############################################################
    fin = readFile(mylist[0])
    recordhash = {}
    dbhash = {}
    index = 0
    for line in fin:
        index += 1
        record = line.strip().split()
        if record[0] not in recordhash:
            recordhash[record[0]] = []
            conn, c = conndbvcf(tmppath, foldernumber, tablename, record[0])
            dbhash[record[0]] = {}
            dbhash[record[0]]["conn"] = conn
            dbhash[record[0]]["c"] = c
        recordhash[record[0]].append(tuple(record))
        if len(recordhash[record[0]]) == nrecord:
            dbhash[record[0]]["c"].executemany("INSERT INTO [{0}] (chromosome, start, Mid, ref, alt) VALUES (?,?,?,?,?)".format(tablename), recordhash[record[0]])
            dbhash[record[0]]["conn"].commit()
            recordhash[record[0]] = []
            # for kk in recordhash[record[0]]:
            #     print(kk)
            #     dbhash[record[0]]["c"].execute("INSERT INTO [{0}] (chromosome, start, Mid, ref, alt) VALUES (?,?,?,?,?)".format(tablename), kk)
            #     dbhash[record[0]]["conn"].commit()
            # recordhash[record[0]] = []
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
    for k in recordhash.keys():
        if len(recordhash[k]) > 0:
            dbhash[k]["c"].executemany("INSERT INTO [{0}] (chromosome, start, Mid, ref, alt) VALUES (?,?,?,?,?)".format(tablename), recordhash[k])
            dbhash[k]["conn"].commit()
        dbhash[k]["c"].close()
        dbhash[k]["conn"].close()
    return index

def conndbbedpe(tmppath, foldernumber, tablename, chromosome):
    tmppath = os.path.join(tmppath, "{0}.{1}.sqlite3".format(foldernumber, chromosome))
    conn = sqlite3.connect(tmppath)
    c = conn.cursor()
    bedpeformat = '''CREATE TABLE IF NOT EXISTS [{0}]
    (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
    chromosome1     CHAR(100)         NOT NULL,
    start1          INT               NOT NULL,
    end1            INT               NOT NULL,
    chromosome2     CHAR(100)         NOT NULL,
    start2          INT               NOT NULL,
    end2            INT               NOT NULL,
    value           FLOAT             NOT NULL,
    barcode         CHAR(100));'''.format(tablename)
    bedpeformatindex = '''CREATE INDEX [idx_{0}] ON [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value, barcode);'''.format(tablename)
    c.execute(bedpeformat)
    c.execute(bedpeformatindex)
    conn.commit()
    return conn, c

def uploadbedpe(mylist, foldernumber, foldpath):
    tmppath = os.path.join(foldpath, "{0}".format(foldernumber))
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)
    print(tmppath)
    ###############################################################
    #                         loop/pet/scpet                      #
    ###############################################################
    fin = readFile(mylist[0])
    recordhash = {}
    dbhash = {}
    index = 0
    barcodedict = {}
    for line in fin:
        record = line.strip().split()
        if record[0] != record[3]:
            continue
        index += 1
        if checkbedpe(record):
            sys.stderr.write("{0} line record is error!, please check\n".format(index))
            sys.stderr.flush()
            sys.exit("{0}\n".format("\t".join(record)))
        # record[6] =  round(float(record[6]), precision)
        if record[0] not in recordhash:
            recordhash[record[0]] = []
            conn, c = conndbbedpe(tmppath, foldernumber, tablename, record[0])
            dbhash[record[0]] = {}
            dbhash[record[0]]["conn"] = conn
            dbhash[record[0]]["c"] = c
        recordhash[record[0]].append(tuple(record))
        if mylist[3] == "mbed" or mylist[3] == "scATAC" or mylist[3] == "scRNA" or mylist[3] == "scPET":
            if record[7] not in barcodedict:
                barcodedict[record[7]] = 0
            barcodedict[record[7]] += 1
        if len(recordhash[record[0]]) == nrecord:
            if mylist[3] == "scPET":
                dbhash[record[0]]["c"].executemany("INSERT INTO [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value, barcode) VALUES (?,?,?,?,?,?,?,?)".format(tablename), recordhash[record[0]])
            else:
                dbhash[record[0]]["c"].executemany("INSERT INTO [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value) VALUES (?,?,?,?,?,?,?)".format(tablename), recordhash[record[0]])
            dbhash[record[0]]["conn"].commit()
            recordhash[record[0]] = []
        if index % 10000000 == 0:
            sys.stdout.write("push 10000000 records ...\n")
            sys.stdout.flush()
    for k in recordhash.keys():
        if len(recordhash[k]) > 0:
            if mylist[3] == "scPET":
                dbhash[k]["c"].executemany("INSERT INTO [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value, barcode) VALUES (?,?,?,?,?,?,?,?)".format(tablename), recordhash[k])
            else:
                dbhash[k]["c"].executemany("INSERT INTO [{0}] (chromosome1, start1, end1, chromosome2, start2, end2, value) VALUES (?,?,?,?,?,?,?)".format(tablename), recordhash[k])
            dbhash[k]["conn"].commit()
        dbhash[k]["c"].close()
        dbhash[k]["conn"].close()

    if mylist[3] == "scPET":
        connttt = sqlite3.connect(os.path.join(tmppath, "{0}.barcodeindex.sqlite3".format(foldernumber)))
        cttt = connttt.cursor()
        cttt.execute('''CREATE TABLE IF NOT EXISTS [{0}]
        (ID            INTEGER PRIMARY KEY   AUTOINCREMENT,
        barcode        CHAR(100)         NOT NULL,
        count          INT               NOT NULL);'''.format(tablename))
        connttt.commit()
        tmplist = []
        barcodelist = list(barcodedict.keys())
        random.shuffle(barcodelist)
        nindex = 0
        for k in barcodelist:
            tmplist.append((k, barcodedict[k]))
            nindex += 1
            if nindex % nrecord == 0:
                cttt.executemany('''INSERT INTO [{0}] (barcode, count) VALUES (?,?)'''.format(tablename), tmplist)
                tmplist = []
        if len(tmplist) > 0:
            cttt.executemany('''INSERT INTO [{0}] (barcode, count) VALUES (?,?)'''.format(tablename), tmplist)
        connttt.commit()
    return index

def upload(args, commands):
    datapath = os.path.join(databasedir, datadir)
    if not os.path.exists(datapath):
        os.makedirs(datapath)
    all_items = os.listdir(datapath)
    folder_names = [int(item) for item in all_items if os.path.isdir(os.path.join(datapath, item))]
    if len(folder_names) == 0:
        curdat1 = 0
    else:
        curdat1 = max(folder_names)
    
    cifdatapath = os.path.join(databasedir, cifdatadir)
    if not os.path.exists(cifdatapath):
        os.makedirs(cifdatapath)
    all_items = os.listdir(cifdatapath)
    folder_names = [int(item) for item in all_items if os.path.isdir(os.path.join(cifdatapath, item))]
    if len(folder_names) == 0:
        curdat2 = 0
    else:
        curdat2 = max(folder_names)

    if curdat1 > curdat2:
        curdat = curdat1
    else:
        curdat = curdat2
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
        deletef       INTEGER               NOT NULL,
        attr          CHAR(255),
        groupID       CHAR(100));'''.format(tablename))
    conn.commit()
    for line in fin:
        if line.startswith("#"):
            fout.write("{0}\tid\n".format(line.strip()))
            continue
        curdat += 1
        tmp = line.strip().split()
        if tmp[3] == "bw" or tmp[3] == "hic" or tmp[3] == "mcool" or tmp[3] == "cif":
            nrecords = uploadbinary(tmp, curdat, datapath)
        elif tmp[3] == "bedGraph" or tmp[3] == "ssbedGraph":
            nrecords = uploadbedgraph(tmp, curdat, datapath)
        elif tmp[3] == "bed6" or tmp[3] == "scRNA" or tmp[3] == "scATAC" or tmp[3] == "mbed" or tmp[3] == "chromhmm" or tmp[3] == "anno":
            nrecords = uploadbed6(tmp, curdat, datapath)
        elif tmp[3] == "vcf":
            nrecords = uploadvcf(tmp, curdat, datapath)
        elif tmp[3] == "bedpe" or tmp[3] == "scPET":
            nrecords = uploadbedpe(tmp, curdat, datapath)
        else:
            sys.exit("file format error")
        c.execute('INSERT INTO [{0}] (sample, format, fileName, assembly, nrecords, dataID, deletef) VALUES (?,?,?,?,?,?,?)'.format(tablename), (tmp[2], tmp[3], tmp[0], tmp[4], nrecords, curdat, 0))
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
    cifdatapath = os.path.join(databasedir, cifdatadir)
    if not os.path.exists(cifdatapath):
        os.makedirs(cifdatapath)
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
    parser_add.add_argument("-f", "--format", required=False, default="AUTO", choices=["AUTO", "bigwig", "hic", "mcool", "cif", "bedgraph", "ssbedgraph", "bed", "vcf", "bedpe", "scRNA", "scATAC", "scPET", "mbed", "anno", "chromhmm"] ,help="Format of input file. If set to AUTO, the program will determine the file format based on the filename suffix. If the inputdirectory parameter is set, it defaults to AUTO [default:AUTO]")
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
