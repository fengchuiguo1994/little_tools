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
import shutil

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
binlength = 50000000

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

def binarytable(tablename):
    return '''CREATE TABLE IF NOT EXISTS `{0}`
            (`ID`         INT           AUTO_INCREMENT, 
            `filePath`    VARCHAR(255)  NOT NULL,
            PRIMARY KEY (ID) )'''.format(tablename)

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

def addTrack(args, commands):
    if args.assembly == spe:
        sys.exit("{0} is key word for assembly information, please change".format(spe))
    if args.loadname == header:
        sys.exit("{0} is key word for table information, please change".format(header))
    try:
        if not os.path.exists(binarypath):
            os.mkdir(binarypath)
        tmpfile = "{0}/tmp{1}tmp".format(binarypath, time.time())
        fout = writeFile(tmpfile)
        fout.write("1")
        fout.close()
        os.remove(tmpfile)
    except (Exception) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while check the database: {0}".format(args.assembly))
        print("Error: {0}. {1}".format(e.message, e.input_value))
        sys.exit("No permission in {0}.\n".format(binarypath))

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
        cursor.execute(binarytable(args.loadname))
        cursor.execute('''INSERT INTO `{0}` (filePath) VALUES ("{1}/{2}")'''.format(args.loadname, binarypath, args.loadname))
        cursor.execute('''INSERT INTO `{0}` (sample,format,fileName,project,nrecords) VALUES ("{1}", "{2}", "{3}", "{4}", "{5}")'''.format(header, args.loadname, args.format, os.path.abspath(args.input), args.project, 1))
        connection.commit()
        cursor.close()
        connection.close()
        shutil.copy(args.input, "{0}/{1}".format(binarypath, args.loadname))
        flag = False
    except (Exception,mysql.connector.Error) as e:
        print("Error: {0}".format(e.__class__.__name__))
        print("Error while creat the table: {0}".format(args.loadname))
        print("Error: {0}. {1}".format(e.message, e.input_value))
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
    if flag:
        sys.exit("please check your parameters")

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

    parser_add = subparsers.add_parser('addTrack',help='add table to database')
    parser_add.add_argument("-i", "--input", required=True, help="the input table")
    parser_add.add_argument("-d", "--assembly", required=True, help="the database/assembly name")
    parser_add.add_argument("-p", "--project", required=True, help = "project name")
    parser_add.add_argument("-n", "--loadname", required=True, help = "the table name")
    parser_add.add_argument('-f', "--format", required=True, choices=["mcool","hic"] ,help="table format")
    parser_add.add_argument("-t", "--threads", default=8, type=int, help = "number of threads [default: 8]")
    parser_add.add_argument("-g", "--groupname", required=False, default=None , help = "the group name, do not use it lightly. May lead to unforeseen errors [default: None]")
    parser_add.set_defaults(func=addTrack)

    commands = sys.argv[1:]
    if ((not commands) or ((commands[0] in ['add', 'remove']) and len(commands) == 1)):
        commands.append('-h')
    args = parser.parse_args(commands)
    return args, commands

if __name__ == "__main__":
    args, commands = getargs()
    args.func(args, commands)