# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/python

import argparse
import sqlite3
import sys
import os
import gzip

# pip install mysql-connector-python
import mysql.connector
from mysql.connector import Error

from multiprocessing import Pool

host = "localhost"
user = "root"
password = "kkltmax8"
spe = "assembly"
precision = 4
header = "header"
nrecord = 5000
binlength = 5000000

class DatabaseInputError(Exception):
    def __init__(self, input_value, message="Invalid input database value"):
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

def upload_data(args):
    data, db_config = args
    db = mysql.connector.connect(**db_config)
    cursor = db.cursor()
    cursor.executemany("INSERT INTO users (name, email) VALUES (%s, %s)", data)
    db.commit()
    cursor.close()
    db.close()
    return len(data)

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

def addGenome(args, commands):
    try:
        # connection = mysql.connector.connect(host = "localhost", user = "root", password = "kkltmax8")
        connection = mysql.connector.connect(host = host, user = user, password = password)
        if connection.is_connected():
            cursor = connection.cursor()
            # 查询数据库
            cursor.execute("SHOW DATABASES")
            databases = [db[0] for db in cursor]
            print(databases)
            if args.assembly in databases:
                raise DatabaseInputError("This database already exists")
            cursor.execute('''CREATE DATABASE IF NOT EXISTS `{0}` DEFAULT CHARSET utf8mb4 COLLATE utf8mb4_unicode_ci'''.format(args.assembly))
            cursor.execute( '''CREATE DATABASE IF NOT EXISTS `{0}` DEFAULT CHARSET utf8mb4 COLLATE utf8mb4_unicode_ci'''.format(spe))
            cursor.execute("use `{0}`".format(spe))
            cursor.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`       INT          AUTO_INCREMENT PRIMARY KEY,
                `assembly`  VARCHAR(50)  NOT NULL) ;'''.format(spe))
            cursor.execute('''INSERT INTO `{0}` (assembly) VALUES ("{1}")'''.format(spe, args.assembly))
            connection.commit()
            print("{0} database created successfully".format(args.assembly))
        connection.rollback()
    except DatabaseInputError as e:
        print("Error while creating database: {0}".format(args.assembly))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        connection.rollback()
    except mysql.connector.Error as err:
        connection.rollback()
    except Error as e:
        print("Error while creating database: ", e)
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
            print("MySQL connection is closed")

def rmGenome(args, commands):
    try:
        # connection = mysql.connector.connect(host = "localhost", user = "root", password = "kkltmax8")
        connection = mysql.connector.connect(host = host, user = user, password = password)
        if connection.is_connected():
            cursor = connection.cursor()
            # 查询数据库
            cursor.execute("SHOW DATABASES")
            databases = [db[0] for db in cursor]
            print(databases)
            if args.assembly not in databases:
                raise DatabaseNotExistsError("This database is not exists")
            cursor.execute('''DROP DATABASE `{0}`'''.format(args.assembly))
            cursor.execute("use `{0}`".format(spe))

            cursor.execute('''SELECT * FROM `{0}` WHERE `{0}` = "{1}";'''.format(spe, args.assembly))
            results = cursor.fetchall()
            for row in results:
                print(row)
                cursor.execute('''DELETE FROM `{0}` WHERE `id` = {1};'''.format(spe,row[0]))
            connection.commit()
            print("{0} database/assembly was deleted.".format(args.assembly))
    except DatabaseNotExistsError as e:
        print("Error while delete database: {0}".format(args.assembly))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        connection.rollback()
    except mysql.connector.Error as err:
        print("please check your parameters")
        connection.rollback()
    except Error as e:
        print("Error while creating database: ", e)
        connection.rollback()
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
            print("MySQL connection is closed")

def insertsize(tablename, inputfile, genome):
    flag = True
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database=genome)
        c = conn.cursor()
        # check table exists
        c.execute("SHOW TABLES")
        tables = [tb[0] for tb in c]
        # print(tables)
        if header in tables:
            c.execute('''SELECT * FROM `{0}` WHERE `format` = "size"'''.format(header))
            results = c.fetchall()
            if len(results) > 0:
                raise SizeDupError("Invalid the size file can only have one")
        
        # create table
        c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
        (`ID`          INT           AUTO_INCREMENT PRIMARY KEY,
        `chromosome`   VARCHAR(50)   NOT NULL,
        `length`       INT           NOT NULL,
        `color`        VARCHAR(10) ) ;'''.format(tablename))
        
        # read input file    
        recordlist = []
        index = 0
        fin = readFile(inputfile)
        for line in fin:
            index += 1
            record = line.strip().split()
            recordlist.append(tuple(record))
            if index % nrecord == 0:
                if len(recordlist[0]) == 3:
                    c.executemany("INSERT INTO `{0}` (chromosome, length, color) VALUES (%s,%s,%s)".format(tablename), recordlist)
                elif len(recordlist[0]) == 2:
                    c.executemany("INSERT INTO `{0}` (chromosome, length) VALUES (%s,%s)".format(tablename), recordlist)
                else:
                    SizeError("Input file formatting error")
                recordlist = []
        if len(recordlist) > 0:
            if len(recordlist[0]) == 3:
                c.executemany("INSERT INTO `{0}` (chromosome, length, color) VALUES (%s,%s,%s)".format(tablename), recordlist)
            elif len(recordlist[0]) == 2:
                c.executemany("INSERT INTO `{0}` (chromosome, length) VALUES (%s,%s)".format(tablename), recordlist)
            else:
                SizeError("Input file formatting error")
        fin.close()
        conn.commit()
    except BedError as e:
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        conn.rollback()
        flag = False
    except mysql.connector.Error as err:
        print("please check your input file")
        conn.rollback()
        flag = False
    except Error as e:
        print("Error load the input file: ", e)
        conn.rollback()
        flag = False
    finally:
        if conn.is_connected():
            c.close()
            conn.close()
            # print("MySQL connection is closed")
    return flag
    
def insertbed6(tablename, inputfile, formatf, genome, chromlist, indexf = None):
    flag = True
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database=genome)
        c = conn.cursor()
        if formatf == "chrband":
            # check table exists
            indexf = None
            c.execute("SHOW TABLES")
            tables = [tb[0] for tb in c]
            # print(tables)
            if header in tables:
                c.execute('''SELECT * FROM `{0}` WHERE `format` = "chrband"'''.format(header))
                results = c.fetchall()
                if len(results) > 0:
                    raise ChrBandDupError("Invalid the chrband file can only have one")
        
        # create table
        if indexf == None:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `chromosome`  VARCHAR(50)   NOT NULL,
            `start`       INT           NOT NULL,
            `end`         INT           NOT NULL,
            `name`        VARCHAR(50)   NOT NULL,
            `score`       FLOAT         NOT NULL,
            `strand`      VARCHAR(1)    NOT NULL,
            `color`       VARCHAR(10),
            `start_bin`   INT           NOT NULL,
            PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED ;'''.format(tablename))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))
        else:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `chromosome`  VARCHAR(50)   NOT NULL,
            `start`       INT           NOT NULL,
            `end`         INT           NOT NULL,
            `name`        VARCHAR(50)   NOT NULL,
            `score`       FLOAT         NOT NULL,
            `strand`      VARCHAR(1)    NOT NULL,
            `color`       VARCHAR(10),
            `start_bin`   INT           NOT NULL,
            PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, indexf))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))

        # read input file    
        recordlist = []
        index = 0
        fin = readFile(inputfile)
        for line in fin:
            record = line.strip().split()
            if record[0] not in chromlist:
                continue
            index += 1
            if checkbed(record):
                raise BedError("This file has error records!")
            record[4] =  round(float(record[4]), precision)
            if record[5] == "-" or record[5] == "0":
                record[5] == "-"
            elif record[5] == "+" or record[5] == "." or record[5] == "1":
                record[5] == "+"
            record.append(int(record[1])//binlength)
            recordlist.append(tuple(record))
            if index % nrecord == 0:
                if len(recordlist[0]) == 8:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, name, score, strand, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                elif len(recordlist[0]) == 7:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, name, score, strand, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                else:
                    BedError("Input file formatting error")
                recordlist = []
        if len(recordlist) > 0:
            if len(recordlist[0]) == 8:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, name, score, strand, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            elif len(recordlist[0]) == 7:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, name, score, strand, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            else:
                BedError("Input file formatting error")
        fin.close()
        conn.commit()
    except BedError as e:
        flag = False
        conn.rollback()
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except mysql.connector.Error as e:
        flag = False
        conn.rollback()
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        print("please check your input file")
    except Error as e:
        flag = False
        conn.rollback()
        print("Error load the input file: ", e)
    finally:
        if flag == False:
            c.execute('''drop table `{0}`'''.format(tablename))
        if conn.is_connected():
            c.close()
            conn.close()
            print("MySQL connection is closed")
    return flag

def insertanno(tablename, inputfile, genome, chromlist, indexf = None):
    flag = True
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database=genome)
        c = conn.cursor()
        # create table
        if indexf == None:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`         INT           AUTO_INCREMENT, 
                `chromosome`  VARCHAR(50)   NOT NULL,
                `start`       INT           NOT NULL,
                `end`         INT           NOT NULL,
                `geneID`      VARCHAR(50)   NOT NULL,
                `geneName`    VARCHAR(50)   NOT NULL,
                `strand`      VARCHAR(1)    NOT NULL,
                `type`        VARCHAR(20)   NOT NULL,
                `color`       VARCHAR(10),
                `start_bin`   INT           NOT NULL,
                PRIMARY KEY (ID, chromosome, start_bin)  ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED ;'''.format(tablename))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))
        else:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`         INT           AUTO_INCREMENT, 
                `chromosome`  VARCHAR(50)   NOT NULL,
                `start`       INT           NOT NULL,
                `end`         INT           NOT NULL,
                `geneID`      VARCHAR(50)   NOT NULL,
                `geneName`    VARCHAR(50)   NOT NULL,
                `strand`      VARCHAR(1)    NOT NULL,
                `type`        VARCHAR(20)   NOT NULL,
                `color`       VARCHAR(10),
                `start_bin`   INT           NOT NULL,
                PRIMARY KEY (ID, chromosome, start_bin)  ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, indexf))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))
        
        # read input file
        recordlist = []
        index = 0
        fin = readFile(inputfile)
        for line in fin:
            record = line.strip().split()
            if record[0] not in chromlist:
                continue
            index += 1
            if checkbed(record):
                raise BedError("This file has error records!")
            if record[5] == "-" or record[5] == "0":
                record[5] == "-"
            elif record[5] == "+" or record[5] == "." or record[5] == "1":
                record[5] == "+"
            record.append(int(record[1])//binlength)
            recordlist.append(tuple(record))
            if index % nrecord == 0:
                if len(recordlist[0]) == 9:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, geneID, geneName, strand, type, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                elif len(recordlist[0]) == 8:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, geneID, geneName, strand, type, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                else:
                    BedError("Input file formatting error")
                recordlist = []
        if len(recordlist) > 0:
            if len(recordlist[0]) == 9:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, geneID, geneName, strand, type, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            elif len(recordlist[0]) == 8:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, geneID, geneName, strand, type, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            else:
                BedError("Input file formatting error")
        fin.close()
        conn.commit()
    except BedError as e:
        flag = False
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except mysql.connector.Error as e:
        flag = False
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        print("please check your input file")
    except Error as e:
        flag = False
        print("Error load the input file: ", e)
    finally:
        if flag == False:
            c.execute('''drop table `{0}`'''.format(tablename))
        if conn.is_connected():
            c.close()
            conn.close()
            print("MySQL connection is closed")
    return flag

def insertbedgraph(tablename, inputfile, genome, chromlist, indexf = None):
    flag = True
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database=genome)
        c = conn.cursor()
        # create table
        if indexf == None:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                    (`ID`         INT           AUTO_INCREMENT,
                    `chromosome`  VARCHAR(50)   NOT NULL,
                    `start`       INT           NOT NULL,
                    `end`         INT           NOT NULL,
                    `value`       FLOAT         NOT NULL,
                    `color`       VARCHAR(10),
                    `start_bin`   INT           NOT NULL,
                    PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ;'''.format(tablename)) # ROW_FORMAT=COMPRESSED
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))
        else:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                    (`ID`         INT           AUTO_INCREMENT,
                    `chromosome`  VARCHAR(50)   NOT NULL,
                    `start`       INT           NOT NULL,
                    `end`         INT           NOT NULL,
                    `value`       FLOAT         NOT NULL,
                    `color`       VARCHAR(10),
                    `start_bin`   INT           NOT NULL,
                    PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB {1}'''.format(tablename, indexf)) # ROW_FORMAT=COMPRESSED
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))
        
        # read input file
        index = 0
        recordlist = []
        fin = readFile(inputfile)
        for line in fin:
            record = line.strip().split()
            if record[0] not in chromlist:
                continue
            index += 1
            if checkbed(record):
                raise BedError("This file has error records!")
            record[3] =  round(float(record[3]), precision)
            record.append(int(record[1])//binlength)
            recordlist.append(tuple(record))
            if index % 1000000 == 0:
                sys.stderr.write("load {0}M records\n".format(index/1000000))
                sys.stderr.flush()
            if index % nrecord == 0:
                if len(recordlist[0]) == 6:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, value, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                elif len(recordlist[0]) == 5:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, value, start_bin) VALUES (%s,%s,%s,%s,%s)".format(tablename), recordlist)
                else:
                    BedError("Input file formatting error")
                recordlist = []
        if len(recordlist) > 0:
            if len(recordlist[0]) == 6:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, value, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            elif len(recordlist[0]) == 5:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, value, start_bin) VALUES (%s,%s,%s,%s,%s)".format(tablename), recordlist)
            else:
                BedError("Input file formatting error")
        fin.close()
        conn.commit()
    except BedError as e:
        flag = False
        conn.rollback()
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except mysql.connector.Error as e:
        flag = False
        conn.rollback()
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        print("please check your input file")
    except Error as e:
        flag = False
        conn.rollback()
        print("Error load the input file: ", e)
    finally:
        if flag == False:
            c.execute('''drop table `{0}`'''.format(tablename))
        if conn.is_connected():
            c.close()
            conn.close()
            print("MySQL connection is closed")
    return flag

def insertssbedgraph(tablename, inputfile, genome, chromlist, indexf = None):
    flag = True
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database=genome)
        c = conn.cursor()
        # create table
        if indexf == None:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`         INT          AUTO_INCREMENT,
                `chromosome`  VARCHAR(50)  NOT NULL,
                `start`       INT          NOT NULL,
                `end`         INT          NOT NULL,
                `value1`      FLOAT        NOT NULL,
                `value2`      FLOAT        NOT NULL,
                `color`       VARCHAR(10),
                `start_bin`   INT          NOT NULL,
                PRIMARY KEY (ID, chromosome, start_bin)  ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED ;'''.format(tablename))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))
        else:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`         INT          AUTO_INCREMENT,
                `chromosome`  VARCHAR(50)  NOT NULL,
                `start`       INT          NOT NULL,
                `end`         INT          NOT NULL,
                `value1`      FLOAT        NOT NULL,
                `value2`      FLOAT        NOT NULL,
                `color`       VARCHAR(10),
                `start_bin`   INT          NOT NULL,
                PRIMARY KEY (ID, chromosome, start_bin)  ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, indexf))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))

        # read input file
        index = 0
        recordlist = []
        fin = readFile(inputfile)
        for line in fin:
            record = line.strip().split()
            if record[0] not in chromlist:
                continue
            index += 1
            if checkbed(record):
                raise BedError("This file has error records!")
            record[3] =  round(float(record[3]), precision)
            record[4] =  round(float(record[4]), precision)
            record.append(int(record[1])//binlength)
            recordlist.append(tuple(record))
            if index % nrecord == 0:
                if len(recordlist[0]) == 7:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, value1, value2, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                elif len(recordlist[0]) == 6:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, value1, value2, start_bin) VALUES (%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                else:
                    BedError("Input file formatting error")
                recordlist = []
        if len(recordlist) > 0:
            if len(recordlist[0]) == 7:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, value1, value2, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            elif len(recordlist[0]) == 6:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, value1, value2, start_bin) VALUES (%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            else:
                BedError("Input file formatting error")
        fin.close()
        conn.commit()
    except BedError as e:
        flag = False
        conn.rollback()
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except mysql.connector.Error as e:
        flag = False
        conn.rollback()
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        print("please check your input file")
    except Error as e:
        flag = False
        conn.rollback()
        print("Error load the input file: ", e)
    finally:
        if flag == False:
            c.execute('''drop table `{0}`'''.format(tablename))
        if conn.is_connected():
            c.close()
            conn.close()
            print("MySQL connection is closed")
    return flag

def insertloop(tablename, inputfile, genome, chromlist, indexf = None):
    flag = True
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database=genome)
        c = conn.cursor()
        # create table
        if indexf == None:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`          INT           AUTO_INCREMENT ,
                `chromosome1`  VARCHAR(50)   NOT NULL,
                `start1`       INT           NOT NULL,
                `end1`         INT           NOT NULL,
                `chromosome2`  VARCHAR(50)   NOT NULL,
                `start2`       INT           NOT NULL,
                `end2`         INT           NOT NULL,
                `value`        FLOAT         NOT NULL,
                `color`        VARCHAR(10),
                `start1_bin`   INT           NOT NULL,
                `start2_bin`   INT           NOT NULL,
                PRIMARY KEY (ID, chromosome1, start1_bin, chromosome2, start2_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED ;'''.format(tablename))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start1, end1, start2, end2); '''.format(tablename))
        else:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`          INT           AUTO_INCREMENT ,
                `chromosome1`  VARCHAR(50)   NOT NULL,
                `start1`       INT           NOT NULL,
                `end1`         INT           NOT NULL,
                `chromosome2`  VARCHAR(50)   NOT NULL,
                `start2`       INT           NOT NULL,
                `end2`         INT           NOT NULL,
                `value`        FLOAT         NOT NULL,
                `color`        VARCHAR(10),
                `start1_bin`   INT           NOT NULL,
                `start2_bin`   INT           NOT NULL,
                PRIMARY KEY (ID, chromosome1, start1_bin, chromosome2, start2_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, indexf))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start1, end1, start2, end2); '''.format(tablename))
        
        # read input file
        index = 0
        recordlist = []
        fin = readFile(inputfile)
        for line in fin:
            record = line.strip().split()
            if record[0] not in chromlist:
                continue
            if record[3] not in chromlist:
                continue
            index += 1
            if checkbed(record):
                raise BedError("This file has error records!")
            if checkbed(record[3:]):
                raise BedError("This file has error records!")
            if checkbedpe(record):
                raise BedpeError("This file has error records!")
            record[6] =  round(float(record[6]), precision)
            record.append(int(record[1])//binlength)
            record.append(int(record[4])//binlength)
            recordlist.append(tuple(record))
            if index % nrecord == 0:
                if len(recordlist[0]) == 10:
                    c.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, value, color, start1_bin, start2_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                elif len(recordlist[0]) == 9:
                    c.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, value, start1_bin, start2_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                else:
                    BedpeError("Input file formatting error")
                recordlist = []
        if len(recordlist) > 0:
            if len(recordlist[0]) == 10:
                c.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, value, color, start1_bin, start2_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            elif len(recordlist[0]) == 9:
                c.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, value, start1_bin, start2_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            else:
                BedpeError("Input file formatting error")
        fin.close()
        conn.commit()
    except BedError as e:
        flag = False
        conn.rollback()
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except BedpeError as e:
        flag = False
        conn.rollback()
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except mysql.connector.Error as e:
        flag = False
        conn.rollback()
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        print("please check your input file")
    except Error as e:
        flag = False
        conn.rollback()
        print("Error load the input file: ", e)
    finally:
        if flag == False:
            c.execute('''drop table `{0}`'''.format(tablename))
        if conn.is_connected():
            c.close()
            conn.close()
            print("MySQL connection is closed")
    return flag

def insertsingleton(tablename, inputfile, genome, chromlist, indexf = None):
    flag = True
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database=genome)
        c = conn.cursor()
        # create table
        if indexf == None:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`         INT          AUTO_INCREMENT,
                `chromosome`  VARCHAR(50)  NOT NULL,
                `start`       INT          NOT NULL,
                `end`         INT          NOT NULL,
                `barcode`     VARCHAR(50)  NOT NULL,
                `value`       INT          NOT NULL,
                `color`       VARCHAR(10),
                `start_bin`   INT          NOT NULL,
                PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED ;'''.format(tablename))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))
        else:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`         INT          AUTO_INCREMENT,
                `chromosome`  VARCHAR(50)  NOT NULL,
                `start`       INT          NOT NULL,
                `end`         INT          NOT NULL,
                `barcode`     VARCHAR(50)  NOT NULL,
                `value`       INT          NOT NULL,
                `color`       VARCHAR(10),
                `start_bin`   INT          NOT NULL,
                PRIMARY KEY (ID, chromosome, start_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, indexf))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome, start, end);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start, end); '''.format(tablename))
        # read input file
        index = 0
        recordlist = []
        fin = readFile(inputfile)
        for line in fin:
            record = line.strip().split()
            if record[0] not in chromlist:
                continue
            index += 1
            if checkbed(record):
                raise BedError("This file has error records!")
            record.append(int(record[1])//binlength)
            recordlist.append(tuple(record))
            if index % 1000000 == 0:
                sys.stderr.write("load {0}M records\n".format(index/1000000))
                sys.stderr.flush()
            if index % nrecord == 0:
                if len(recordlist[0]) == 7:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, barcode, value, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                elif len(recordlist[0]) == 6:
                    c.executemany("INSERT INTO `{0}` (chromosome, start, end, barcode, value, start_bin) VALUES (%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                else:
                    BedError("Input file formatting error")
                recordlist = []
        if len(recordlist) > 0:
            if len(recordlist[0]) == 7:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, barcode, value, color, start_bin) VALUES (%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            elif len(recordlist[0]) == 6:
                c.executemany("INSERT INTO `{0}` (chromosome, start, end, barcode, value, start_bin) VALUES (%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            else:
                BedError("Input file formatting error")
        fin.close()
        conn.commit()
    except BedError as e:
        flag = False
        conn.rollback()
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except mysql.connector.Error as e:
        flag = False
        conn.rollback()
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        print("please check your input file")
    except Error as e:
        flag = False
        conn.rollback()
        print("Error load the input file: ", e)
    finally:
        if flag == False:
            c.execute('''drop table `{0}`'''.format(tablename))
        if conn.is_connected():
            c.close()
            conn.close()
            print("MySQL connection is closed")
    return flag

def insertscPET(tablename, inputfile, genome, chromlist, indexf = None):
    flag = True
    try:
        conn = mysql.connector.connect(host = host, user = user, password = password, database=genome)
        c = conn.cursor()
        # create table
        if indexf == None:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`          INT           AUTO_INCREMENT,
                `chromosome1`  VARCHAR(50)   NOT NULL,
                `start1`       INT           NOT NULL,
                `end1`         INT           NOT NULL,
                `chromosome2`  VARCHAR(50)   NOT NULL,
                `start2`       INT           NOT NULL,
                `end2`         INT           NOT NULL,
                `barcode`      VARCHAR(50)   NOT NULL,
                `value`        INT           NOT NULL,
                `color`        VARCHAR(10),
                `start1_bin`   INT           NOT NULL,
                `start2_bin`   INT           NOT NULL,
                PRIMARY KEY (ID, chromosome1, start1_bin, chromosome2, start2_bin)  ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED ;'''.format(tablename))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start1, end1, start2, end2); '''.format(tablename))
        else:
            c.execute('''CREATE TABLE IF NOT EXISTS `{0}`
                (`ID`          INT           AUTO_INCREMENT,
                `chromosome1`  VARCHAR(50)   NOT NULL,
                `start1`       INT           NOT NULL,
                `end1`         INT           NOT NULL,
                `chromosome2`  VARCHAR(50)   NOT NULL,
                `start2`       INT           NOT NULL,
                `end2`         INT           NOT NULL,
                `barcode`      VARCHAR(50)   NOT NULL,
                `value`        INT           NOT NULL,
                `color`        VARCHAR(10),
                `start1_bin`   INT           NOT NULL,
                `start2_bin`   INT           NOT NULL,
                PRIMARY KEY (ID, chromosome1, start1_bin, chromosome2, start2_bin) ) ENGINE=InnoDB ROW_FORMAT=COMPRESSED {1}'''.format(tablename, indexf))
            c.execute('''CREATE INDEX `idx_{0}` ON `{0}` (chromosome1, start1, end1, chromosome2, start2, end2);'''.format(tablename))
            # c.execute('''ALTER TABLE `{0}` ADD SPATIAL INDEX `idx_{0}_spatial` (start1, end1, start2, end2); '''.format(tablename))
        
        # read input file
        index = 0
        recordlist = []
        fin = readFile(inputfile)
        for line in fin:
            record = line.strip().split()
            if record[0] not in chromlist:
                continue
            if record[3] not in chromlist:
                continue
            index += 1
            if checkbed(record):
                raise BedError("This file has error records!")
            if checkbed(record[3:]):
                raise BedError("This file has error records!")
            if checkbedpe(record):
                raise BedpeError("This file has error records!")
            record.append(int(record[1])//binlength)
            record.append(int(record[4])//binlength)
            recordlist.append(tuple(record))
            if index % nrecord == 0:
                if len(recordlist[0]) == 11:
                    c.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, barcode, value, color, start1_bin, start2_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                elif len(recordlist[0]) == 10:
                    c.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, barcode, value, start1_bin, start2_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
                else:
                    BedpeError("Input file formatting error")
                recordlist = []
        if len(recordlist) > 0:
            if len(recordlist[0]) == 11:
                c.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, barcode, value, color, start1_bin, start2_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            elif len(recordlist[0]) == 10:
                c.executemany("INSERT INTO `{0}` (chromosome1, start1, end1, chromosome2, start2, end2, barcode, value, start1_bin, start2_bin) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(tablename), recordlist)
            else:
                BedpeError("Input file formatting error")
        fin.close()
        conn.commit()
    except BedError as e:
        flag = False
        conn.rollback()
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except BedpeError as e:
        flag = False
        conn.rollback()
        print("Error while read the input file: {0}".format(inputfile))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
    except mysql.connector.Error as e:
        flag = False
        conn.rollback()
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        print("please check your input file")
    except Error as e:
        flag = False
        conn.rollback()
        print("Error load the input file: ", e)
    finally:
        if flag == False:
            c.execute('''drop table `{0}`'''.format(tablename))
        if conn.is_connected():
            c.close()
            conn.close()
            print("MySQL connection is closed")
    return flag
          
def addTrack(args, commands):
    flag = False
    if args.loadname == header:
        sys.exit("The {0} is a summary table of all the information and should not be used by the user, so please use another name.\n".format(header))
    try:
        # connection = mysql.connector.connect(host = "localhost", user = "root", password = "kkltmax8")
        connection = mysql.connector.connect(host = host, user = user, password = password)
        if connection.is_connected():
            cursor = connection.cursor()
            # 查询数据库
            cursor.execute("SHOW DATABASES")
            databases = [db[0] for db in cursor]
            print(databases)
            if args.assembly not in databases:
                raise DatabaseNotExistsError("This database is not exists")
            cursor.execute('''USE `{0}`'''.format(args.assembly))
            cursor.execute("SHOW TABLES")
            tables = [tb[0] for tb in cursor]
            print(tables)
            if args.loadname in tables:
                raise TableInputError("This table already exists")
    except DatabaseNotExistsError as err:
        connection.rollback()
        print("Error: {0}".format(err))
        flag = True
    except TableInputError as err:
        connection.rollback()
        print("Error: {0}".format(err))
        flag = True
    except mysql.connector.Error as err:
        # 发生错误时回滚事务
        connection.rollback()
        print("Error: {0}".format(err))
        flag = True
    finally:
        # 关闭游标和连接
        cursor.close()
        connection.close()
    if flag:
        sys.exit("please check your parameters\n")

    # push data
    flag = False
    # size file (genome size file, samtools faidx)
    if args.format == "size":
        flag = insertsize(args.loadname, args.input, args.assembly)
    else:
        connection = mysql.connector.connect(host = host, user = user, password = password, database=args.assembly)
        cursor = connection.cursor()
        cursor.execute("SHOW tables;")
        tmptables = [tbe[0] for tbe in cursor]
        if header not in tmptables:
            cursor.close()
            connection.close()
            sys.exit("Must have size table1\n")
        cursor.execute('''SELECT sample FROM `{0}` WHERE format = "size";'''.format(header))
        tmptables = [tbe[0] for tbe in cursor]
        print(tmptables)
        if len(tmptables) == 0:
            cursor.close()
            connection.close()
            sys.exit("Must have size table2\n")
        elif len(tmptables) > 1:
            cursor.close()
            connection.close()
            sys.exit("there are two size table\n")
        else:
            cursor.execute('''SELECT * FROM `{0}`;'''.format(tmptables[0]))
            cslist = [cs for cs in cursor] # chromosome list
            chromlist = set()
            partbed = "PARTITION BY LIST COLUMNS(chromosome, start_bin) ("
            partbedpe = "PARTITION BY LIST COLUMNS(chromosome1, start1_bin, chromosome2, start2_bin) ("
            for index1,(id1,chromosome1,length1,color1) in enumerate(cslist):
                chromlist.add(chromosome1)
                nbin1 = 0
                while nbin1 * binlength <= length1:  # start: nbin*binlength, end: (nbin+1)*binlength
                    partbed += ''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1)
                    # print(''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1))
                    for index2,(id2,chromosome2,length2,color2) in enumerate(cslist):
                        nbin2 = 0
                        while nbin2 * binlength <= length2:  # start: nbin*binlength, end: (nbin+1)*binlength
                            partbedpe += ''' PARTITION p{0}_{1}_{2}_{3} VALUES IN (("{4}", {5}, "{6}", {7})), '''.format(index1, nbin1, index2, nbin2, chromosome1, nbin1, chromosome2, nbin2)
                            nbin2 += 1
                            # print(''' PARTITION p{0}_{1}_{2}_{3} VALUES IN (("{4}", {5}, "{6}", {7})), '''.format(index1, nbin1, index2, nbin2, chromosome1, nbin1, chromosome2, nbin2))
                    nbin1 += 1    
            partbed = "{0} );".format(partbed[:-2])
            partbedpe = "{0} );".format(partbedpe[:-2])
            print("generate split table info")
        
        # chrband file #
        if args.format == "chrband":
            flag = insertbed6(args.loadname, args.input, args.format, args.assembly, partbed)
        # bed6/chromhmm file #
        elif args.format == "bed6" or args.format == "chromhmm":
            flag = insertbed6(args.loadname, args.input, args.format, args.assembly, partbed)
        # anno
        elif args.format == "anno":
            flag = insertanno(args.loadname, args.input, args.assembly, partbed)
        # bedgraph (bulk RNA-Seq/ATAC-Seq...)
        elif args.format == "bedgraph" or args.format == "RNAbedgraph" or args.format == "ATACbedgraph":
            flag = insertbedgraph(args.loadname, args.input, args.assembly, partbed)
        # ssbedgraph (bulk ssRNA-Seq)
        elif args.format == "ssbedgraph":
            flag = insertssbedgraph(args.loadname, args.input, args.assembly, partbed)
        # loop (ChIA-PET cluster/Hi-C loop)    
        elif args.format == "loop":
            flag = insertloop(args.loadname, args.input, args.assembly, partbedpe)
        # singleton (scRNA/scATAC/scBS) multi (multi-way interaction/ChIA-Drop/SPRITE/Pore-C
        elif args.format == "scRNA" or args.format == "scATAC" or args.format == "multi" or args.format == "singleton":
            flag = insertsingleton(args.loadname, args.input, args.assembly, partbed)
        # scPET (scHi-C/scPET/scChIATAC) 
        elif args.format == "scPET":
            flag = insertscPET(args.loadname, args.input, args.assembly, partbedpe)
        else:
            print(args.format)
            sys.exit("error format\n")

    # add information into header table          
    if flag:
        try:
            connection = mysql.connector.connect(host = host, user = user, password = password, database = args.assembly)
            cursor = connection.cursor()
            cursor.execute('''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`        INT            AUTO_INCREMENT  PRIMARY KEY,
            `sample`     VARCHAR(50)    NOT NULL,
            `format`     VARCHAR(20)    NOT NULL,
            `fileName`   VARCHAR(255)   NOT NULL,
            `project`    VARCHAR(50)    NOT NULL,
            `group`      VARCHAR(50) );'''.format(header))
            if args.groupname != None:
                cursor.execute('''INSERT INTO `{0}` (sample,format,fileName,project,group) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}")'''.format(header, args.loadname, args.format, os.path.abspath(args.input), args.project, args.groupname))
            else:
                cursor.execute('''INSERT INTO `{0}` (sample,format,fileName,project) VALUES ("{1}", "{2}", "{3}", "{4}")'''.format(header, args.loadname, args.format, os.path.abspath(args.input), args.project))
            connection.commit()
        except mysql.connector.Error as err:
            print("please check your input file", err)
            connection.rollback()
        except Error as e:
            print("Error load the input file: ", e)
            connection.rollback()
        finally:
            if connection.is_connected():
                cursor.close()
                connection.close()
                # print("MySQL connection is closed")
    print("Store data success")

def rmTrack(args, commands):
    if args.loadname == header:
        sys.exit("The {0} is a summary table of all the information and should not be used by the user, so please use another name.\n".format(header))
    try:
        # connection = mysql.connector.connect(host = "localhost", user = "root", password = "kkltmax8")
        connection = mysql.connector.connect(host = host, user = user, password = password)
        if connection.is_connected():
            cursor = connection.cursor()
            # 查询数据库
            cursor.execute("SHOW DATABASES")
            databases = [db[0] for db in cursor]
            print(databases)
            if args.assembly not in databases:
                raise DatabaseNotExistsError("This database is not exists")
            cursor.execute('''USE `{0}`'''.format(args.assembly))
            cursor.execute("SHOW TABLES")
            tables = [tb[0] for tb in cursor]
            print(tables)
            if args.loadname not in tables:
                raise TableInputError("This table is not exists")
        
        cursor.execute('''DROP table `{0}`'''.format(args.loadname))
        cursor.execute('''SELECT * FROM `{0}` WHERE sample = "{1}";'''.format(header, args.loadname))
        results = cursor.fetchall()
        for row in results:
            print("delete table {0}".format(row))
            cursor.execute('''DELETE FROM `{0}` WHERE `id` = {1};'''.format(header,row[0]))
        connection.commit()
    except DatabaseNotExistsError as e:
        print("Error while delete table: {0} in database: {1}".format(args.loadname, args.assembly))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        connection.rollback()
    except TableInputError as e:
        print("Error while delete table: {0}".format(args.loadname))
        print("Error: {0}. The input value was: {1}".format(e.message, e.input_value))
        connection.rollback()
    except mysql.connector.Error as err:
        print("please check your parameters")
        connection.rollback()
    except Error as e:
        print("Error while creating database: ", e)
        connection.rollback()
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
            print("MySQL connection is closed")
    print("{0} table was deleted.".format(args.loadname))
        
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

    flag = False
    try:
        # connection = mysql.connector.connect(host = "localhost", user = "root", password = "kkltmax8")
        connection = mysql.connector.connect(host = host, user = user, password = password)
        if connection.is_connected():
            cursor = connection.cursor()
            # 查询数据库
            cursor.execute("SHOW DATABASES")
            databases = [db[0] for db in cursor]
            print(databases)
            if args.assembly not in databases:
                raise DatabaseNotExistsError("This database is not exists")
            cursor.execute('''USE `{0}`'''.format(args.assembly))
            cursor.execute("SHOW TABLES")
            tables = [tb[0] for tb in cursor]
            print(tables)
            if header in tables:
                # print('''SELECT * FROM `{0}` WHERE `group` = "{1}" AND `sample` = "{2}";'''.format(header, args.groupname, args.loadname))
                cursor.execute('''SELECT * FROM `{0}` WHERE `group` = "{1}" AND `sample` = "{2}";'''.format(header, args.groupname, args.loadname))
                tmp = [tb for tb in cursor]
                if len(tmp) > 0:
                    raise TableInputError("This table already exists. {0}".format("\n".join(tmp)))
    except DatabaseNotExistsError as err:
        connection.rollback()
        print("Error: {0}".format(err))
        flag = True
    except TableInputError as err:
        connection.rollback()
        print("Error: {0}".format(err))
        flag = True
    except mysql.connector.Error as err:
        connection.rollback()
        print("Error: {0}".format(err))
        flag = True
    finally:
        # 关闭游标和连接
        cursor.close()
        connection.close()
    if flag:
        sys.exit("please check your parameters\n")

    # split table
    connection = mysql.connector.connect(host = host, user = user, password = password, database=args.assembly)
    cursor = connection.cursor()
    cursor.execute("SHOW tables;")
    tmptables = [tbe[0] for tbe in cursor]
    if header not in tmptables:
        cursor.close()
        connection.close()
        sys.exit("Must have size table3\n")
    cursor.execute('''SELECT sample FROM `{0}` WHERE format = "size";'''.format(header))
    tmptables = [tbe[0] for tbe in cursor]
    if len(tmptables) == 0:
        cursor.close()
        connection.close()
        sys.exit("Must have size table4\n")
    elif len(tmptables) > 1:
        cursor.close()
        connection.close()
        sys.exit("there are two size table\n")
    else:
        cursor.execute('''SELECT * FROM `{0}`;'''.format(tmptables[0]))
        cslist = [cs for cs in cursor] # chromosome list
        chromlist = set()
        partbed = "PARTITION BY LIST COLUMNS(chromosome, start_bin) ("
        partbedpe = "PARTITION BY LIST COLUMNS(chromosome1, start1_bin, chromosome2, start2_bin) ("
        for index1,(id1,chromosome1,length1,color1) in enumerate(cslist):
            chromlist.add(chromosome1)
            nbin1 = 0
            while nbin1 * binlength <= length1:  # start: nbin*binlength, end: (nbin+1)*binlength
                partbed += ''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1)
                # print(''' PARTITION p{0}_{1} VALUES IN (("{2}",{3})), '''.format(index1, nbin1, chromosome1, nbin1))
                for index2,(id2,chromosome2,length2,color2) in enumerate(cslist):
                    nbin2 = 0
                    while nbin2 * binlength <= length2:  # start: nbin*binlength, end: (nbin+1)*binlength
                        partbedpe += ''' PARTITION p{0}_{1}_{2}_{3} VALUES IN (("{4}", {5}, "{6}", {7})), '''.format(index1, nbin1, index2, nbin2, chromosome1, nbin1, chromosome2, nbin2)
                        nbin2 += 1
                        # print(''' PARTITION p{0}_{1}_{2}_{3} VALUES IN (("{4}", {5}, "{6}", {7})), '''.format(index1, nbin1, index2, nbin2, chromosome1, nbin1, chromosome2, nbin2))
                nbin1 += 1    
        partbed = "{0} );".format(partbed[:-2])
        partbedpe = "{0} );".format(partbedpe[:-2])
        print("generate split table info")

    # push data
    flags = []
    for iii,jjj in zip(inputs,formats):
        if jjj == "RNAbedgraph":
            flag = insertbedgraph("{0}.{1}.RNAbedgraph".format(args.loadname, args.groupname), iii, args.assembly, partbed)
        elif  jjj == "ATACbedgraph":
            flag = insertbedgraph("{0}.{1}.ATACbedgraph".format(args.loadname, args.groupname), iii, args.assembly, partbed)
        elif jjj == "loop":
            flag = insertloop("{0}.{1}.loop".format(args.loadname, args.groupname), iii, args.assembly, partbedpe)
        elif jjj == "scRNA":
            flag = insertsingleton("{0}.{1}.scRNA".format(args.loadname, args.groupname), iii, args.assembly, partbed)
        elif jjj == "scATAC":
            flag = insertsingleton("{0}.{1}.scATAC".format(args.loadname, args.groupname), iii, args.assembly, partbed)
        elif jjj == "scPET":
            flag = insertscPET("{0}.{1}.scPET".format(args.loadname, args.groupname), iii, args.assembly, partbedpe)
        else:
            print(iii,jjj)
            sys.exit("error format\n")
        flags.append(flag)

    # add information into header table          
    if False not in flags:
        try:
            connection = mysql.connector.connect(host = host, user = user, password = password, database = args.assembly)
            cursor = connection.cursor()
            cursor.execute('''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`        INT            AUTO_INCREMENT  PRIMARY KEY,
            `sample`     VARCHAR(50)    NOT NULL,
            `format`     VARCHAR(20)    NOT NULL,
            `fileName`   VARCHAR(255)   NOT NULL,
            `project`    VARCHAR(50)    NOT NULL,
            `group`      VARCHAR(50) );'''.format(header))
            for iii,jjj in zip(inputs,formats):
                cursor.execute('''INSERT INTO `{0}` (`sample`,`format`,`fileName`,`project`,`group`) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}")'''.format(header, "{0}.{1}.{2}".format(args.loadname, args.groupname, jjj), jjj, os.path.abspath(iii), args.project, args.groupname))
            connection.commit()
        except mysql.connector.Error as err:
            print("please check your input file", err)
            connection.rollback()
        except Error as e:
            print("Error load the input file: ", e)
            connection.rollback()
        finally:
            if connection.is_connected():
                cursor.close()
                connection.close()
                print("MySQL connection is closed")
    print("Store data success")

def plotAssembly(args, commands):
    connection = mysql.connector.connect(host = host, user = user, password = password)
    if connection.is_connected():
        cursor = connection.cursor()
        # 查询数据库
        cursor.execute("SHOW DATABASES")
        databases = [db[0] for db in cursor]
        if args.assembly not in databases:
            sys.exit("This assembly/database is not exists\n")
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
    parser_add.set_defaults(func=addGenome)

    parser_add = subparsers.add_parser('rmGenome',help='remove database/assembly')
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.set_defaults(func=rmGenome)

    parser_add = subparsers.add_parser('addTrack',help='add table to database')
    parser_add.add_argument("-i", "--input", required=True, help="the input table")
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.add_argument("-p", "--project", required=True, help = "project name")
    parser_add.add_argument("-n", "--loadname", required=True, help = "the table name")
    parser_add.add_argument("-g", "--groupname", required=False, default=None , help = "the group name, do not use it lightly. May lead to unforeseen errors [default: None]")
    parser_add.add_argument('-f', "--format", default="bedgraph" ,choices=["bedgraph","ssbedgraph","loop","scRNA","scATAC","scPET","multi","anno","chromhmm","size","chrband","bed6"] ,help="table format [default:bedgraph]")
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
    parser_add.add_argument('-f', "--format", required=True ,help="the input table format [format1,file2,file3,file4(,file5,file6)] [scRNA,scATAC,scPET,RNAbedgraph,ATACbedgraph,loop]")
    parser_add.set_defaults(func=addMultiTrack)

    parser_add = subparsers.add_parser('plotAssembly',help='plot assembly/databases information')
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.set_defaults(func=plotAssembly)

    commands = sys.argv[1:]
    if ((not commands) or ((commands[0] in ['add', 'remove']) and len(commands) == 1)):
        commands.append('-h')
    args = parser.parse_args(commands)
    return args, commands

if __name__ == "__main__":
    args, commands = getargs()
    args.func(args, commands)