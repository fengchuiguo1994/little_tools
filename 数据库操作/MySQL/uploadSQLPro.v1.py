# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys
import os
import gzip
import argparse
import logging
import multiprocessing as mp
# from contextlib import contextmanager
# from typing import Dict
import random
from functools import partial
import time

# pip install mysql-connector-python
from mysql.connector import MySQLConnection
from mysql.connector.cursor import MySQLCursor
import mysql.connector

logging.basicConfig(level=logging.INFO,format = '%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

host = "localhost"
user = "root"
password = "kkltmax8"
autocommit = False

spe = "assembly"
precision = 4
header = "header"
nrecord = 5000
binlength = 5000000

binarypath = "/home/binarypath"

bedset = set(["chrband", "bed6", "chromhmm", "anno", "bedgraph", "RNAbedgraph", "ATACbedgraph", "ssbedgraph", "scRNA", "scATAC", "multi", "singleton"])
bedsetbed6 = set(["chrband", "bed6", "chromhmm", "scRNA", "scATAC", "multi", "singleton"]) # chr s e name score strand color bin
bedsetbed8 = set(["anno"]) # chr s e name score strand genename type color bin
bedsetbedbedgraph = set(["bedgraph", "RNAbedgraph", "ATACbedgraph"]) # chr s e value color bin
bedsetbedssbedgraph = set(["ssbedgraph"]) # # chr s e value1 value2 color1 color2 bin
bedpeset = set(["loop", "scPET"])

def headertable(header):
    return '''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`        INT            AUTO_INCREMENT  PRIMARY KEY,
            `sample`     VARCHAR(50)    NOT NULL,
            `format`     VARCHAR(20)    NOT NULL,
            `fileName`   VARCHAR(255)   NOT NULL,
            `project`    VARCHAR(50)    NOT NULL,
            `nrecords`   INT            NOT NULL,
            `group`      VARCHAR(50) );'''.format(header)

def bed6table(tablename, partbed):
    return '''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `chromosome`  VARCHAR(50)   NOT NULL,
            `start`       INT           NOT NULL,
            `end`         INT           NOT NULL,
            `name`        VARCHAR(50)   NOT NULL,
            `score`       FLOAT         NOT NULL,
            `strand`      VARCHAR(1)    NOT NULL,
            `color`       VARCHAR(10)   NOT NULL,
            `start_bin`   INT           NOT NULL,
            PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, partbed)

def bed8table(tablename, partbed):
    return '''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `chromosome`  VARCHAR(50)   NOT NULL,
            `start`       INT           NOT NULL,
            `end`         INT           NOT NULL,
            `name`        VARCHAR(50)   NOT NULL,
            `score`       FLOAT         NOT NULL,
            `strand`      VARCHAR(1)    NOT NULL,
            `geneName`    VARCHAR(50)   NOT NULL,
            `type`        VARCHAR(20)   NOT NULL,
            `color`       VARCHAR(10)   NOT NULL,
            `start_bin`   INT           NOT NULL,
            PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, partbed)

def bedgraphtable(tablename, partbed):
    return '''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `chromosome`  VARCHAR(50)   NOT NULL,
            `start`       INT           NOT NULL,
            `end`         INT           NOT NULL,
            `value`       FLOAT         NOT NULL,
            `color`       VARCHAR(10)   NOT NULL,
            `start_bin`   INT           NOT NULL,
            PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, partbed)

def ssbedgraphtable(tablename, partbed):
    return '''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `chromosome`  VARCHAR(50)   NOT NULL,
            `start`       INT           NOT NULL,
            `end`         INT           NOT NULL,
            `value1`      FLOAT         NOT NULL,
            `value2`      FLOAT         NOT NULL,
            `color1`      VARCHAR(10)   NOT NULL,
            `color2`      VARCHAR(10)   NOT NULL,
            `start_bin`   INT           NOT NULL,
            PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, partbed)

def looptable(tablename, partbed):
    return '''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `chromosome1`  VARCHAR(50)   NOT NULL,
            `start1`       INT           NOT NULL,
            `end1`         INT           NOT NULL,
            `chromosome2`  VARCHAR(50)   NOT NULL,
            `start2`       INT           NOT NULL,
            `end2`         INT           NOT NULL,
            `value`        FLOAT         NOT NULL,
            `color`        VARCHAR(10)   NOT NULL,
            `start1_bin`   INT           NOT NULL,
            PRIMARY KEY (ID, chromosome1, start1_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, partbed)

def scpettable(tablename, partbed):
    return '''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `chromosome1`  VARCHAR(50)   NOT NULL,
            `start1`       INT           NOT NULL,
            `end1`         INT           NOT NULL,
            `chromosome2`  VARCHAR(50)   NOT NULL,
            `start2`       INT           NOT NULL,
            `end2`         INT           NOT NULL,
            `barcode`      VARCHAR(50)   NOT NULL,
            `value`        FLOAT         NOT NULL,
            `color`        VARCHAR(10)   NOT NULL,
            `start1_bin`   INT           NOT NULL,
            PRIMARY KEY (ID, chromosome1, start1_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, partbed)

class DatabaseInputError(Exception):
    def __init__(self, input_value, message="Invalid input database value"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)

class KeyWordInputError(Exception):
    def __init__(self, input_value, message="Invalid input value"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)

class DatabaseNotExistsError(Exception):
    def __init__(self, input_value, message="Database is not exists"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)

class TableInputError(Exception):
    def __init__(self, input_value, message="Invalid input table value"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)

class BedError(Exception):
    def __init__(self, input_value, message="Invalid Start >= end"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)

class BedpeError(Exception):
    def __init__(self, input_value, message="Invalid must anchor1 < anchor2"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)

class SizeError(Exception):
    def __init__(self, input_value, message="Invalid file must have two column: chromosome, length"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)
    
class SizeDupError(Exception):
    def __init__(self, input_value, message="Invalid the size file can only have one"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)

class ChrBandDupError(Exception):
    def __init__(self, input_value, message="Invalid the chrband file can only have one"):
        self.input_value = input_value
        self.message = message
        super().__init__(self.message)

def checkbed(mylist):
    if mylist[1].isdigit() and mylist[2].isdigit():
        if int(mylist[1]) < int(mylist[2]):
            return False
    return True

def is_oppo_number(s):
    try:
        aa = float(s)
        if aa >= 0:
            return True
        else:
            return False
    except ValueError:
        return False
    
def checksize(mylist):
    if len(mylist) == 2 or len(mylist) == 3:
        if mylist[1].isdigit():
            return False
    return True

def compareChar(a1, a2):
    if a1 > a2:
        return 1
    elif a1 == a2:
        return 0
    else:
        return -1
    
def checkbedpe(mylist):
    if mylist[1].isdigit() and mylist[2].isdigit() and mylist[4].isdigit() and mylist[5].isdigit():
        if int(mylist[1]) < int(mylist[2]) and int(mylist[4]) < int(mylist[5]):    
            result = compareChar(mylist[0], mylist[3])
            if result == 1:
                return True
            elif result == -1:
                return False
            else:
                # if int(mylist[1]) < int(mylist[4]):
                if int(mylist[1]) <= int(mylist[4]): # unsolve bug
                    return False
    return True

def readFile(infile):
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin

def addGenome(args, commands):
    try:
        connection = mysql.connector.connect(host = host, user = user, password = password, autocommit = True)
        cursor = connection.cursor()
        ######## check
        if args.assembly == spe:
            raise KeyWordInputError("{0} is key word for assembly information, please change".format(spe))
        if "{0}.size".format(args.assembly) == header:
            raise KeyWordInputError("{0} is key word for table information, please change".format(header))
        cursor.execute('''show databases;''')
        databases = [db[0] for db in cursor]    
        if args.assembly in databases:
            raise DatabaseInputError("The {0} database already exists".format(args.assembly))

        ########## add genome
        cursor.execute('''CREATE DATABASE `{0}` DEFAULT CHARSET utf8mb4 COLLATE utf8mb4_unicode_ci'''.format(args.assembly))
        cursor.execute( '''CREATE DATABASE IF NOT EXISTS `{0}` DEFAULT CHARSET utf8mb4 COLLATE utf8mb4_unicode_ci'''.format(spe))
        cursor.execute("use `{0}`".format(spe))
        cursor.execute('''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`       INT          AUTO_INCREMENT PRIMARY KEY,
            `assembly`  VARCHAR(50)  NOT NULL) ;'''.format(spe))
        cursor.execute('''INSERT INTO `{0}` (assembly) VALUES ("{1}")'''.format(spe, args.assembly))

        cursor.execute("use `{0}`".format(args.assembly))
        cursor.execute('''CREATE TABLE `{0}`
        (`ID`          INT           AUTO_INCREMENT PRIMARY KEY,
        `chromosome`   VARCHAR(50)   NOT NULL,
        `length`       INT           NOT NULL) ;'''.format("{0}.size".format(args.assembly)))
        cursor.execute(headertable(header))
        connection.commit()

        #########  add size
        recordlist = []
        index = 0
        fin = readFile(args.sizefile)
        for line in fin:
            index += 1
            record = line.strip().split()
            if checksize(record):
                raise SizeError("Input file formatting error")
            recordlist.append(tuple(record))
            if index % nrecord == 0:
                cursor.executemany("INSERT INTO `{0}.size` (chromosome, length) VALUES (%s,%s)".format(args.assembly), recordlist)
                recordlist = []
        if len(recordlist) > 0:
            cursor.executemany("INSERT INTO `{0}.size` (chromosome, length) VALUES (%s,%s)".format(args.assembly), recordlist)
        fin.close()
        cursor.execute('''INSERT INTO `{0}` (sample,format,fileName,project,nrecords) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}")'''.format(header, "{0}.size".format(args.assembly), "size", os.path.abspath(args.sizefile), "anno", index))
        connection.commit()
        cursor.close()
        connection.close()

    except (KeyWordInputError,SizeError,mysql.connector.Error,DatabaseInputError,Exception) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while creating database: {0}".format(args.assembly))
        print("Error: {0}. {1}".format(e.message, e.input_value))
        if e.__class__.__name__ == "KeyWordInputError" or e.__class__.__name__ == "DatabaseInputError":
            pass
        else:
            try:
                cursor.execute("use `{0}`".format(spe))
                cursor.execute('''DELETE FROM `{0}` where assembly = "{1}";'''.format(spe, args.assembly))
                connection.commit()
            except Exception:
                pass
            try:
                cursor.execute("use `{0}`".format(spe))
                cursor.execute('''DROP DATABASE `{0}`'''.format(args.assembly))
                connection.commit()
            except Exception:
                pass
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
    
def rmGenome(args, commands):
    try:
        connection = mysql.connector.connect(host = host, user = user, password = password)
        cursor = connection.cursor()
        cursor.execute("SHOW DATABASES")
        databases = [db[0] for db in cursor]
        if args.assembly not in databases:
            raise DatabaseNotExistsError("This database is not exists")
        cursor.execute('''DROP DATABASE `{0}`'''.format(args.assembly))
        cursor.execute("use `{0}`".format(spe))
        cursor.execute('''DELETE FROM `{0}` WHERE `{0}` = "{1}";'''.format(spe, args.assembly))
        connection.commit()
        logger.info("{0} database/assembly was deleted.".format(args.assembly))
        cursor.close()
        connection.close()
    except (DatabaseNotExistsError,mysql.connector.Error,Exception) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while delete database: {0}".format(args.assembly))
        print("Error: {0}. {1}".format(e.message, e.input_value))
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()

def addTrack(args, commands):
    if args.assembly == spe:
        sys.exit("{0} is key word for assembly information, please change".format(spe))
    if args.loadname == header:
        sys.exit("{0} is key word for table information, please change".format(header))
    flag = True
    try:
        connection = mysql.connector.connect(host = host, user = user, password = password)
        cursor = connection.cursor()
        ########## check
        cursor.execute('''show databases;''')
        databases = [db[0] for db in cursor]    
        if args.assembly not in databases:
            raise DatabaseNotExistsError("The {0} database is not exists".format(args.assembly))
        cursor.execute("use `{0}`".format(spe))
        cursor.execute('''SELECT `{0}` FROM `{0}` WHERE `{0}` = "{1}";'''.format(spe, args.assembly))
        databases = [db[0] for db in cursor]
        if args.assembly not in databases:
            raise DatabaseNotExistsError("The {0} database is not exists".format(args.assembly))
        
        cursor.execute('''USE `{0}`'''.format(args.assembly))
        cursor.execute("SHOW TABLES")
        tables = [tb[0] for tb in cursor]
        if args.loadname in tables:
            raise TableInputError("The {0} table already exists".format(args.loadname))
        cursor.execute('''SELECT sample FROM `{0}` WHERE sample = "{1}";'''.format(header, args.loadname))
        tables = [tb[0] for tb in cursor]
        if args.loadname in tables:
            raise TableInputError("The {0} table already exists".format(args.loadname))

        if args.format == "chrband":
            cursor.execute('''SELECT sample FROM `{0}` WHERE format = "chrband";'''.format(header))
            tmptables = [tbe[0] for tbe in cursor]
            if len(tmptables) > 0:
                raise ChrBandDupError("The [{0}] table already exists".format("/".join(tmptables)))

        cursor.execute('''SELECT sample FROM `{0}` WHERE format = "size";'''.format(header))
        tmptables = [tbe[0] for tbe in cursor]
        cursor.execute('''SELECT * FROM `{0}`;'''.format(tmptables[0]))
        cslist = [cs for cs in cursor] # chromosome list
        chromlist = set()
        partbed = "PARTITION BY LIST COLUMNS(chromosome, start_bin) ("
        partbed2 = "PARTITION BY LIST COLUMNS(chromosome1, start1_bin) ("
        # partbedpe = "PARTITION BY LIST COLUMNS(chromosome1, start1_bin, chromosome2, start2_bin) ("
        for index1,(id1,chromosome1,length1) in enumerate(cslist):
            chromlist.add(chromosome1)
            nbin1 = 0
            while nbin1 * binlength <= length1:  # start: nbin*binlength, end: (nbin+1)*binlength
                partbed += ''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1)
                partbed2 += ''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1)
                # print(''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1))
                # for index2,(id2,chromosome2,length2) in enumerate(cslist):
                #     nbin2 = 0
                #     while nbin2 * binlength <= length2:  # start: nbin*binlength, end: (nbin+1)*binlength
                #         partbedpe += ''' PARTITION p{0}_{1}_{2}_{3} VALUES IN (("{4}", {5}, "{6}", {7})), '''.format(index1, nbin1, index2, nbin2, chromosome1, nbin1, chromosome2, nbin2)
                #         nbin2 += 1
                #         # print(''' PARTITION p{0}_{1}_{2}_{3} VALUES IN (("{4}", {5}, "{6}", {7})), '''.format(index1, nbin1, index2, nbin2, chromosome1, nbin1, chromosome2, nbin2))
                nbin1 += 1    
        partbed = "{0} );".format(partbed[:-2])
        partbed2 = "{0} );".format(partbed2[:-2])
        # partbedpe = "{0} );".format(partbedpe[:-2])
        cursor.close()
        connection.close()
        flag = False
    except (DatabaseNotExistsError,TableInputError,mysql.connector.Error,Exception) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while check the database: {0}".format(args.assembly))
        print("Error: {0}. {1}".format(e.message, e.input_value))
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
    if flag:
        sys.exit("please check your parameters")
    
    flag = True
    try:
        connection = mysql.connector.connect(host = host, user = user, password = password, database = args.assembly)
        cursor = connection.cursor()
        ############ load data
        if args.format in bedsetbed6:
            cursor.execute(bed6table(args.loadname, partbed))
            cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(args.loadname))
            # cursor.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(args.loadname))
        elif args.format in bedsetbed8:
            cursor.execute(bed8table(args.loadname, partbed))
            cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(args.loadname))
            # cursor.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(args.loadname))
        elif args.format in bedsetbedbedgraph:
            cursor.execute(bedgraphtable(args.loadname, partbed))
            cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(args.loadname))
            # cursor.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(args.loadname))
        elif args.format in bedsetbedssbedgraph:
            cursor.execute(ssbedgraphtable(args.loadname, partbed))
            cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(args.loadname))
            # cursor.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(args.loadname))
        elif args.format == "loop":
            cursor.execute(looptable(args.loadname, partbed2))
            cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(args.loadname))
            # cursor.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start1, end1, start2, end2); '''.format(args.loadname))
        elif args.format == "scPET":
            cursor.execute(scpettable(args.loadname, partbed2))
            cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(args.loadname))
            # cursor.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start1, end1, start2, end2); '''.format(args.loadname))
        else:
            sys.exit("miss1\n")

        totalrecord = {"total":0}
        if args.format == "scPET" or args.format == "scRNA" or args.format == "scATAC":
            barcodedict = {}
            chunk_gen = get_chunks(nrecord, args.input, args.format, chromlist, totalrecord, barcodedict)
        else:
            chunk_gen = get_chunks(nrecord, args.input, args.format, chromlist, totalrecord)
        nthread = args.threads
        if nthread > mp.cpu_count()-1:
            nthread = mp.cpu_count()-1
        pool = mp.Pool(nthread)
        paraworker = partial(workerMulti, args.assembly, args.format, args.loadname)
        results = pool.imap(paraworker, chunk_gen)
        pool.close()
        pool.join()
        
        ############# index 
        if args.groupname == None:
            cursor.execute('''INSERT INTO `{0}` (sample,format,fileName,project,nrecords) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}")'''.format(header, args.loadname, args.format, os.path.abspath(args.input), args.project, totalrecord["total"]))
        else:
            cursor.execute('''INSERT INTO `{0}` (sample,format,fileName,project,nrecords,group) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}", "{6}")'''.format(header, args.loadname, args.format, os.path.abspath(args.input), args.project, totalrecord["total"], args.groupname))
        connection.commit()

        if args.format == "scPET" or args.format == "scRNA" or args.format == "scATAC":
            cursor.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`       INT          AUTO_INCREMENT PRIMARY KEY,
                `barcode`   VARCHAR(50)  NOT NULL,
                `count`     INT          NOT NULL) ;'''.format("{0}.barcodeindex".format(args.loadname)))
            index = 0
            recordlist = []
            barcodelist = list(barcodedict.keys())
            random.shuffle(barcodelist)
            for k in barcodelist:
                index += 1
                recordlist.append((k, barcodedict[k]))
                if index % nrecord == 0:
                    cursor.executemany('''INSERT INTO `{0}` (barcode, count) VALUES (%s,%s)'''.format("{0}.barcodeindex".format(args.loadname)), recordlist)
                    recordlist = []
            if len(recordlist) > 0:
                cursor.executemany('''INSERT INTO `{0}` (barcode, count) VALUES (%s,%s)'''.format("{0}.barcodeindex".format(args.loadname)), recordlist)
            connection.commit()

            if args.groupname == None:
                cursor.execute('''INSERT INTO `{0}` (sample,format,fileName,project,nrecords) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}")'''.format(header, "{0}.barcodeindex".format(args.loadname), "{0}index".format(args.format), os.path.abspath(args.input), args.project, index))
            else:
                cursor.execute('''INSERT INTO `{0}` (sample,format,fileName,project,nrecords,group) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}", "{6}")'''.format(header, "{0}.barcodeindex".format(args.loadname), "{0}index".format(args.format), os.path.abspath(args.input), args.project, index, args.groupname))
            connection.commit()

        connection.commit()
        cursor.close()
        connection.close()
        flag = False
    except (BedError,BedpeError,Exception,mysql.connector.Error) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while creat the table: {0}".format(args.loadname))
        print("Error: {0}. {1}".format(e.message, e.input_value))
        try:
            cursor.execute(''' drop table `{0}`'''.format(args.loadname))
            connection.commit()
        except Exception:
            pass
        try:
            cursor.execute(''' drop table `{0}`'''.format("{0}.barcodeindex".format(args.loadname)))
            connection.commit()
        except Exception:
            pass
        try:
            cursor.execute(''' delete from `{0}` where sample = "{1}"'''.format(header, args.loadname))
            connection.commit()
        except Exception:
            pass
        try:
            cursor.execute(''' delete from `{0}` where sample = "{1}"'''.format(header, "{0}.barcodeindex".format(args.loadname)))
            connection.commit()
        except Exception:
            pass
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
    if flag:
        sys.exit("please check your parameters")

def get_chunks(batch_size, infile, informat, chromset, totalrecord, barcodedict = None):
    try:
        begintime = time.time()
        currenttime = begintime
        recordlist = []
        index = 0
        fin = readFile(infile)
        for line in fin:
            record = line.strip().split()
            if informat in bedset:
                if record[0] not in chromset:
                    continue
                if checkbed(record):
                    sys.stderr.write("Error line {0}: {1}\n".format(index,record))
                    sys.stderr.flush()
                    raise BedError("This file has error records!")
                if informat == "scRNA" or informat == "scATAC":
                    if len(record) == 5:
                        record.append("+")
                if informat in bedsetbed6 or informat in bedsetbed8:
                    # record[4] =  round(float(record[4]), precision)
                    if record[5] == "-" or record[5] == "0":
                        record[5] == "-"
                    elif record[5] == "+" or record[5] == "." or record[5] == "1":
                        record[5] == "+"
                    if len(record) == 6:
                        record.append("-")
                elif informat in bedsetbedbedgraph:
                    if len(record) == 4:
                        record.append("-")
                elif informat == "ssbedgraph":
                    if len(record) == 5:
                        record.append("-")
                        record.append("-")
                record.append(int(record[1])//binlength)
                recordlist.append(tuple(record))
                if barcodedict != None:
                    if record[3] not in barcodedict:
                        barcodedict[record[3]] = 0
                    barcodedict[record[3]] += 1
            elif informat in bedpeset:
                if record[0] not in chromset or record[3] not in chromset:
                    continue
                if checkbedpe(record):
                    sys.stderr.write("Error line {0}: {1}\n".format(index,record))
                    sys.stderr.flush()
                    raise BedpeError("This file has error records!")
                if informat == "loop":
                    if len(record) == 7:
                        record.append("-")
                if informat == "scPET":
                    if len(record) == 7:
                        record.append(1)
                    if len(record) == 8:
                        record.append("-")
                record.append(int(record[1])//binlength)
                recordlist.append(tuple(record))
                if barcodedict != None:
                    if record[6] not in barcodedict:
                        barcodedict[record[6]] = 0
                    barcodedict[record[6]] += 1
            else:
                sys.exit("your format not allow")
            index += 1
            if index % 5000000 == 0:
                print("5M records in {0} has been loaded, use time: {1:.3f} sec.".format(infile, time.time() - currenttime))
                currenttime = time.time()
            if index % batch_size == 0:
                yield recordlist
                recordlist = []
        if recordlist:
            yield recordlist
        totalrecord["total"] = index
        print("{0} records in {1} has been loaded, use time: {2:.3f} sec.".format(index, infile ,time.time() - begintime))
    except (BedError,BedpeError,Exception) as e:
        raise e
    finally:
        fin.close()

def worker(args, chunk_gen):
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database = args.assembly)
        cursor = conn.cursor()
        if args.format in bedsetbed6:
            cursor.executemany("INSERT INTO `{0}` (chromosome, start, end, name, score, strand, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)".format(args.loadname), chunk_gen)
        elif args.format in bedsetbed8:
            cursor.executemany("INSERT INTO `{0}` (chromosome, start, end, name, score, strand, geneName, type, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(args.loadname), chunk_gen)
        elif args.format in bedsetbedbedgraph:
            cursor.executemany("INSERT INTO `{0}` (chromosome, start, end, value, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s)".format(args.loadname), chunk_gen)
        elif args.format in bedsetbedssbedgraph:
            cursor.executemany("INSERT INTO `{0}` (chromosome, start, end, value1, value2, color1, color2, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)".format(args.loadname), chunk_gen)
        elif args.format == "loop":
            cursor.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, value, color, start1_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(args.loadname), chunk_gen)
        elif args.format == "scPET":
            cursor.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, barcode, value, color, start1_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(args.loadname), chunk_gen)
        else:
            sys.exit("miss6")
        conn.commit()
        cursor.close()
        conn.close()
    except mysql.connector.Error as e:
        print(e)
        raise e
    finally:
        if conn.is_connected():
            cursor.close()
            conn.close()

def workerMulti(assemblyname, formatname, uploadname, chunk_gen):
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database = assemblyname)
        cursor = conn.cursor()
        if formatname in bedsetbed6:
            cursor.executemany("INSERT INTO `{0}` (chromosome, start, end, name, score, strand, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)".format(uploadname), chunk_gen)
        elif formatname in bedsetbed8:
            cursor.executemany("INSERT INTO `{0}` (chromosome, start, end, name, score, strand, geneName, type, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(uploadname), chunk_gen)
        elif formatname in bedsetbedbedgraph:
            cursor.executemany("INSERT INTO `{0}` (chromosome, start, end, value, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s)".format(uploadname), chunk_gen)
        elif formatname in bedsetbedssbedgraph:
            cursor.executemany("INSERT INTO `{0}` (chromosome, start, end, value1, value2, color1, color2, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)".format(uploadname), chunk_gen)
        elif formatname == "loop":
            cursor.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, value, color, start1_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(uploadname), chunk_gen)
        elif formatname == "scPET":
            cursor.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, barcode, value, color, start1_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(uploadname), chunk_gen)
        else:
            sys.exit("miss7")
        conn.commit()
        cursor.close()
        conn.close()
    except (mysql.connector.Error, Exception) as e:
        raise e
    finally:
        if conn.is_connected():
            cursor.close()
            conn.close()

def rmTrack(args, commands):
    if args.assembly == spe:
        sys.exit("{0} is key word for assembly information, please change".format(spe))
    if args.loadname == header:
        sys.exit("{0} is key word for table information, please change".format(header))
    flag = True
    try:
        connection = mysql.connector.connect(host = host, user = user, password = password)
        cursor = connection.cursor()
        cursor.execute("SHOW DATABASES")
        databases = [db[0] for db in cursor]
        print(databases)
        if args.assembly not in databases:
            raise DatabaseNotExistsError("This database is not exists")
        cursor.execute('''USE `{0}`'''.format(args.assembly))
        cursor.execute("SHOW TABLES")
        tables = [tb[0] for tb in cursor]
        if args.loadname not in tables:
            raise TableInputError("This table is not exists")
        
        cursor.execute('''SELECT * FROM `{0}` WHERE sample = "{1}";'''.format(header, args.loadname))
        results = cursor.fetchall()
        for row in results:
            print("delete table {0}".format(row))
            cursor.execute('''DELETE FROM `{0}` WHERE `id` = {1};'''.format(header,row[0]))
            if row[2] in ["scRNA","scATAC","scPET"]:
                cursor.execute('''DELETE FROM `{0}` WHERE `sample` = "{1}.barcodeindex";'''.format(header, row[1]))
                cursor.execute('''DROP table `{0}.barcodeindex`'''.format(args.loadname))
        cursor.execute('''DROP table `{0}`'''.format(args.loadname))
        connection.commit()
        cursor.close()
        connection.close()
        print("{0} table was deleted.".format(args.loadname))
    except (DatabaseNotExistsError,TableInputError,mysql.connector.Error) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while creating database: {0}".format(args.assembly))
        print("Error: {0}. {1}".format(e.message, e.input_value))
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
        
def addMultiTrack(args, commands):
    allowFormats = set(["scRNA","scATAC","scPET","RNAbedgraph","ATACbedgraph","loop"])
    inputs = args.input.split(",")
    formats = args.format.split(",")
    if len(inputs) != len(formats):
        sys.exit("must have the same number of files and formats\n")
    if len(set(formats)) != len(formats):
        sys.exit("limited to 1 per format\n")
    flag = False
    for ii in formats:
        if ii not in allowFormats:
            flag = True
            break
    if flag:
        sys.exit("the format must be chosen from these [scRNA,scATAC,scPET,RNAbedgraph,ATACbedgraph,loop] \n")

    ############### check
    flag = True
    try:
        connection = mysql.connector.connect(host = host, user = user, password = password)
        cursor = connection.cursor()
        cursor.execute("SHOW DATABASES")
        databases = [db[0] for db in cursor]
        if args.assembly not in databases:
            raise DatabaseNotExistsError("This database is not exists")
        cursor.execute('''USE `{0}`'''.format(args.assembly))
        cursor.execute("SHOW TABLES")
        tables = [tb[0] for tb in cursor]
        if header in tables:
            cursor.execute('''SELECT * FROM `{0}` WHERE `group` = "{1}" AND `sample` = "{2}";'''.format(header, args.groupname, args.loadname))
            tmp = [tb for tb in cursor]
            if len(tmp) > 0:
                raise TableInputError("This table already exists. {0}".format("\n".join(tmp)))

        cursor.execute('''SELECT sample FROM `{0}` WHERE format = "size";'''.format(header))
        tmptables = [tbe[0] for tbe in cursor]
        cursor.execute('''SELECT * FROM `{0}`;'''.format(tmptables[0]))
        cslist = [cs for cs in cursor] # chromosome list
        chromlist = set()
        partbed = "PARTITION BY LIST COLUMNS(chromosome, start_bin) ("
        partbed2 = "PARTITION BY LIST COLUMNS(chromosome1, start1_bin) ("
        # partbedpe = "PARTITION BY LIST COLUMNS(chromosome1, start1_bin, chromosome2, start2_bin) ("
        for index1,(id1,chromosome1,length1) in enumerate(cslist):
            chromlist.add(chromosome1)
            nbin1 = 0
            while nbin1 * binlength <= length1:  # start: nbin*binlength, end: (nbin+1)*binlength
                partbed += ''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1)
                partbed2 += ''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1)
                # print(''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1))
                # for index2,(id2,chromosome2,length2) in enumerate(cslist):
                #     nbin2 = 0
                #     while nbin2 * binlength <= length2:  # start: nbin*binlength, end: (nbin+1)*binlength
                #         partbedpe += ''' PARTITION p{0}_{1}_{2}_{3} VALUES IN (("{4}", {5}, "{6}", {7})), '''.format(index1, nbin1, index2, nbin2, chromosome1, nbin1, chromosome2, nbin2)
                #         nbin2 += 1
                #         # print(''' PARTITION p{0}_{1}_{2}_{3} VALUES IN (("{4}", {5}, "{6}", {7})), '''.format(index1, nbin1, index2, nbin2, chromosome1, nbin1, chromosome2, nbin2))
                nbin1 += 1    
        partbed = "{0} );".format(partbed[:-2])
        partbed2 = "{0} );".format(partbed2[:-2])
        # partbedpe = "{0} );".format(partbedpe[:-2])
        cursor.close()
        connection.close()
        flag = False
    except (DatabaseNotExistsError,TableInputError,mysql.connector.Error,Exception) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while check database: {0}".format(args.assembly))
        print("Error: {0}. {1}".format(e.message, e.input_value))
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
    if flag:
        sys.exit("please check your parameters\n")

    ############ load data
    flag = set()
    try:
        connection = mysql.connector.connect(host = host, user = user, password = password, database = args.assembly)
        cursor = connection.cursor()
        totalrecordall = {"RNAbedgraph":0, "ATACbedgraph":0, "loop":0, "scRNA":0, "scATAC":0, "scPET":0}
        barcodedictall = {"scRNA":None, "scATAC":None, "scPET":None}
        combinebarcode = {}
        combinen = 0
        for iii,jjj in zip(inputs,formats):
            totalrecord = {"total":0}
            barcodedict = {}
            loadnamegroupname = "{0}.{1}".format(args.loadname, args.groupname)
            if jjj == "RNAbedgraph":
                loadname = "{0}.{1}.RNAbedgraph".format(args.loadname, args.groupname)
                cursor.execute(bedgraphtable(loadname, partbed))
                cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(loadname))
                chunk_gen = get_chunks(nrecord, iii, jjj, chromlist, totalrecord)
            elif  jjj == "ATACbedgraph":
                loadname = "{0}.{1}.ATACbedgraph".format(args.loadname, args.groupname)
                cursor.execute(bedgraphtable(loadname, partbed))
                cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(loadname))
                chunk_gen = get_chunks(nrecord, iii, jjj, chromlist, totalrecord)
            elif jjj == "loop":
                loadname = "{0}.{1}.loop".format(args.loadname, args.groupname)
                cursor.execute(looptable(loadname, partbed2))
                cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(loadname))
                chunk_gen = get_chunks(nrecord, iii, jjj, chromlist, totalrecord)
            elif jjj == "scRNA":
                loadname = "{0}.{1}.scRNA".format(args.loadname, args.groupname)
                cursor.execute(bed6table(loadname, partbed))
                cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(loadname))
                chunk_gen = get_chunks(nrecord, iii, jjj, chromlist, totalrecord, barcodedict)
            elif jjj == "scATAC":
                loadname = "{0}.{1}.scATAC".format(args.loadname, args.groupname)
                cursor.execute(bed6table(loadname, partbed))
                cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(loadname))
                chunk_gen = get_chunks(nrecord, iii, jjj, chromlist, totalrecord, barcodedict)
            elif jjj == "scPET":
                loadname = "{0}.{1}.scPET".format(args.loadname, args.groupname)
                cursor.execute(scpettable(loadname, partbed2))
                cursor.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(loadname))
                chunk_gen = get_chunks(nrecord, iii, jjj, chromlist, totalrecord, barcodedict)
            else:
                print(iii,jjj)
                sys.exit("error format\n")

            nthread = args.threads
            if nthread > mp.cpu_count()-1:
                nthread = mp.cpu_count()-1
            pool = mp.Pool(nthread)
            paraworker = partial(workerMulti, args.assembly, jjj, loadname)
            results = pool.imap(paraworker, chunk_gen)
            pool.close()
            pool.join()

            cursor.execute('''INSERT INTO `{0}` (`sample`,`format`,`fileName`,`project`,`nrecords`,`group`) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}", "{6}")'''.format(header, loadname, jjj, os.path.abspath(iii), args.project, totalrecord["total"], args.groupname))
            connection.commit()

            ############## index 
            if jjj == "RNAbedgraph":
                totalrecordall["RNAbedgraph"] = totalrecord["total"]
            elif  jjj == "ATACbedgraph":
                totalrecordall["ATACbedgraph"] = totalrecord["total"]
            elif jjj == "loop":
                totalrecordall["loop"] = totalrecord["total"]
            elif jjj == "scRNA":
                totalrecordall["scRNA"] = totalrecord["total"]
                barcodedictall["scRNA"] = barcodedict
            elif jjj == "scATAC":
                totalrecordall["scATAC"] = totalrecord["total"]
                barcodedictall["scATAC"] = barcodedict
            elif jjj == "scPET":
                totalrecordall["scPET"] = totalrecord["total"]
                barcodedictall["scPET"] = barcodedict
            else:
                print(iii,jjj)
                sys.exit("error format\n")
            
            if jjj == "scPET" or jjj == "scRNA" or jjj == "scATAC":
                combinen += 1
                cursor.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                    (`ID`       INT          AUTO_INCREMENT PRIMARY KEY,
                    `barcode`   VARCHAR(50)  NOT NULL,
                    `count`     INT          NOT NULL) ;'''.format("{0}.barcodeindex".format(loadname)))
                index = 0
                recordlist = []
                barcodelist = list(barcodedict.keys())
                random.shuffle(barcodelist)
                for k in barcodelist:
                    if k not in combinebarcode:
                        combinebarcode[k] = 0
                    combinebarcode[k] += 1
                    index += 1
                    recordlist.append((k, barcodedict[k]))
                    if index % nrecord == 0:
                        cursor.executemany('''INSERT INTO `{0}` (barcode, count) VALUES (%s,%s)'''.format("{0}.barcodeindex".format(loadname)), recordlist)
                        recordlist = []
                if len(recordlist) > 0:
                    cursor.executemany('''INSERT INTO `{0}` (barcode, count) VALUES (%s,%s)'''.format("{0}.barcodeindex".format(loadname)), recordlist)
                connection.commit()

                cursor.execute('''INSERT INTO `{0}` (`sample`,`format`,`fileName`,`project`,`nrecords`,`group`) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}", "{6}")'''.format(header, "{0}.barcodeindex".format(loadname), "{0}index".format(jjj), os.path.abspath(iii), args.project, index, args.groupname))
                connection.commit()
            
        if combinen > 0:
            cursor.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`       INT          AUTO_INCREMENT PRIMARY KEY,
                `barcode`   VARCHAR(50)  NOT NULL,
                `count`     INT          NOT NULL) ;'''.format("{0}.joint.barcodeindex".format(loadnamegroupname)))
            barcodelist = []
            for k in combinebarcode:
                if combinebarcode[k] == combinen:
                    barcodelist.append(k)
            index = 0
            recordlist = []
            random.shuffle(barcodelist)
            for k in barcodelist:
                index += 1
                recordlist.append((k, 0))
                if index % nrecord == 0:
                    cursor.executemany('''INSERT INTO `{0}` (barcode, count) VALUES (%s,%s)'''.format("{0}.joint.barcodeindex".format(loadnamegroupname)), recordlist)
                    recordlist = []
            if len(recordlist) > 0:
                cursor.executemany('''INSERT INTO `{0}` (barcode, count) VALUES (%s,%s)'''.format("{0}.joint.barcodeindex".format(loadnamegroupname)), recordlist)
            connection.commit()

            cursor.execute('''INSERT INTO `{0}` (`sample`,`format`,`fileName`,`project`,`nrecords`,`group`) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}", "{6}")'''.format(header, "{0}.joint.barcodeindex".format(loadnamegroupname), "jointindex".format(jjj), "no", args.project, index, args.groupname))
            connection.commit()
        cursor.close()
        connection.close()
    except (BedError,BedpeError,Exception,mysql.connector.Error) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while creat the table: {0}.{1}".format(args.loadname, args.groupname))
        print("Error: {0}. {1}".format(e.message, e.input_value))
        for k in ["RNAbedgraph","ATACbedgraph","loop","scRNA","scATAC","scPET"]:    
            try:
                cursor.execute(''' drop table `{0}.{1}.{2}`'''.format(args.loadname, args.groupname,k))
                connection.commit()
            except Exception:
                pass
            try:
                cursor.execute(''' drop table `{0}.{1}.{2}.barcodeindex`'''.format(args.loadname, args.groupname,k))
                connection.commit()
            except Exception:
                pass
            try:
                cursor.execute(''' delete from `{0}` where sample = "{1}.{2}.{3}"'''.format(header, args.loadname, args.groupname,k))
                connection.commit()
            except Exception:
                pass
            try:
                cursor.execute(''' delete from `{0}` where sample = "{1}.{2}.{3}.barcodeindex"'''.format(header, args.loadname, args.groupname,k))
                connection.commit()
            except Exception:
                pass
        try:
            cursor.execute(''' drop table `{0}.{1}.joint.barcodeindex`'''.format(args.loadname, args.groupname))
            connection.commit()
        except Exception:
            pass
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()

def plotGenome(args, commands):
    connection = mysql.connector.connect(host = host, user = user, password = password)
    cursor = connection.cursor()
    if args.assembly == None:
        cursor.execute("SHOW DATABASES")
        databases = [db[0] for db in cursor]
        if spe not in databases:
            print("There are no databases")
        else:
            cursor.execute('''USE `{0}`;'''.format(spe))
            cursor.execute("SHOW TABLES;")
            tables = [tb[0] for tb in cursor]
            if spe not in tables:
                print("There are no databases")
            else:
                cursor.execute('''SELECT assembly FROM `{0}`;'''.format(spe))
                databases = [db[0] for db in cursor]
                print("There are {0} databases;".format("/".join(databases)))
    else:
        cursor.execute("SHOW DATABASES")
        databases = [db[0] for db in cursor]
        if args.assembly not in databases:
            sys.exit("This assembly/database is not exists\n")
        if spe not in databases:
            sys.exit("This assembly/database is not exists\n")
        else:
            cursor.execute('''USE `{0}`;'''.format(spe))
            cursor.execute("SHOW TABLES;")
            tables = [tb[0] for tb in cursor]
            if spe not in tables:
                sys.exit("This assembly/database is not exists\n")
            else:
                cursor.execute('''SELECT assembly FROM `{0}`;'''.format(spe))
                databases = [db[0] for db in cursor]
                if args.assembly in databases:
                    print("database: {0}".format(args.assembly))
                    cursor.execute('''USE `{0}`;'''.format(args.assembly))
                    cursor.execute("SHOW TABLES;")
                    tables = [tb[0] for tb in cursor]
                    if header not in tables:
                        sys.exit("There are no table\n")
                    cursor.execute('''SELECT * FROM `{0}`;'''.format(header))
                    tmptables = [tbe for tbe in cursor]
                    # columns = [column for column in cursor.description]
                    columns = [column[0] for column in cursor.description]
                    if len(tmptables) == 0:
                        print("There are no table")
                    else:
                        print("There are {0} tables: ".format(len(tmptables)))
                        print("\t".join(columns))
                        for ii in tmptables:
                            print('\t'.join(str(item) for item in ii))

def getargs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Program: sqlControllers
Version: 1.0.0
Usage: sqlControllers <command> [options]
Commands:
    addGenome  add the database/assemble.
    add        add table to database (if it does not contain the database/table, it is created).
    remove     remove table.
==========================================================================


""")
    parser.add_argument('--version', action='version',version='%(prog)s {0} by xyhuang'.format('1.0.0'))
    subparsers = parser.add_subparsers(title='subcommands',description='valid subcommands',help='') 
    
    ## add subcommand 
    parser_add = subparsers.add_parser('addGenome',help='add database/assembly')
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.add_argument("-s", "--sizefile", required=True, help="the genome size file")
    parser_add.set_defaults(func=addGenome)

    parser_add = subparsers.add_parser('rmGenome',help='remove database/assembly')
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.set_defaults(func=rmGenome)

    parser_add = subparsers.add_parser('addTrack',help='add table to database')
    parser_add.add_argument("-i", "--input", required=True, help="the input table")
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.add_argument("-p", "--project", required=True, help = "project name")
    parser_add.add_argument("-n", "--loadname", required=True, help = "the table name")
    parser_add.add_argument('-f', "--format", required=True, choices=["bedgraph","ssbedgraph","loop","scRNA","scATAC","scPET","multi","anno","chromhmm","chrband","bed6","RNAbedgraph","ATACbedgraph"] ,help="table format")
    parser_add.add_argument("-t", "--threads", default=8, type=int, help = "number of threads [default: 8]")
    parser_add.add_argument("-g", "--groupname", required=False, default=None , help = "the group name, do not use it lightly. May lead to unforeseen errors [default: None]")
    parser_add.set_defaults(func=addTrack)

    parser_add = subparsers.add_parser('rmTrack',help='add table to database')
    parser_add.add_argument("-n", "--loadname", required=True, help = "the table name")
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.add_argument("-p", "--project", required=True, help = "project name")
    parser_add.set_defaults(func=rmTrack)

    parser_add = subparsers.add_parser('addMultiTrack',help='add tri-omics table to database')
    parser_add.add_argument("-i", "--input", required=True, help="the input table [file1,file2,file3,file4(,file5,file6)]")
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.add_argument("-p", "--project", required=True, help = "project name")
    parser_add.add_argument("-n", "--loadname", required=True, help = "the table name [loadname.groupname.scRNA,loadname.groupname.scATAC,loadname.groupname.scPET,loadname.groupname.RNAbedgraph,loadname.groupname.ATACbedgraph,loadname.groupname.loop]")
    parser_add.add_argument("-g", "--groupname", required=True, help = "the group name")
    parser_add.add_argument("-t", "--threads", default=8, type=int, help = "number of threads [default: 8]")
    parser_add.add_argument('-f', "--format", required=True ,help="the input table format [format1,file2,file3,file4(,file5,file6)] [scRNA,scATAC,scPET,RNAbedgraph,ATACbedgraph,loop]")
    parser_add.set_defaults(func=addMultiTrack)

    parser_add = subparsers.add_parser('plotGenome',help='plot assembly/databases information')
    parser_add.add_argument("-d", "--assembly", default=None, help="the database/assembly name. if not set, it will print all database [default: None]")
    parser_add.set_defaults(func=plotGenome)

    commands = sys.argv[1:]
    if ((not commands) or ((commands[0] in ['add', 'remove']) and len(commands) == 1)):
        commands.append('-h')
    args = parser.parse_args(commands)
    return args, commands

if __name__ == "__main__":
    args, commands = getargs()
    args.func(args, commands)
