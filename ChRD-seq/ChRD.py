import argparse
import sys
import time
import os
import subprocess
import re
import functools
import gzip

def cmp(a, b):
    if b < a:
        return 1
    if a < b:
        return -1
    return 0

def mycmp(o1,o2):
    if isinstance(o1,str) and isinstance(o2,str): 
        return cmp(o1,o2)
    elif isinstance(o1,int) and isinstance(o2,int): 
        return cmp(o1,o2)
    elif isinstance(o1,str) and isinstance(o2,int): 
        return 1
    else:
        return -1

def flash(args,commands):
    flog = open(args.logfile,'a')
    flog.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " use flash combine the R1 and R2 start ...\n")
    flog.flush()
    be = time.time()
    returncode,returnresult = subprocess.getstatusoutput("flash --threads {0} -M {1} --output-prefix {2} --output-directory {3} {4} {5}".format(args.threads,args.maxlength,args.prefix,args.outdir,args.inputR1,args.inputR2))
    if returncode != 0:
        sys.stderr.write("[ERROR]: failed to run flash : \n{0}\n".format(returnresult))
        sys.stderr.flush()
    for i in returnresult.split("\n"):
        if re.search(r"Total pairs:\s*(\d+)",i):
            tol = int(re.search(r"Total pairs:\s*(\d+)",i)[1])
            flog.write("Total pairs:{0}\n".format(tol))
        if re.search(r"Combined pairs:\s*(\d+)",i):
            com = int(re.search(r"Combined pairs:\s*(\d+)",i)[1])
            flog.write("Combined pairs:{0}\n".format(com))
        if re.search(r"Uncombined pairs:\s*(\d+)",i):
            uncom = int(re.search(r"Uncombined pairs:\s*(\d+)",i)[1])
            flog.write("Uncombined pairs:{0}\n".format(uncom))
    flog.write("Combined pairs percentage:{0:.5f}\n".format(com/tol))
    flog.write("Uncombined pairs percentage:{0:.5f}\n".format(uncom/tol))
    flog.write("result file: combine file:{0}/{1}.extendedFrags.fastq, not combine file1:{0}/{1}.notCombined_1.fastq, not combine file2:{0}/{1}.notCombined_2.fastq, hist file:{0}/{1}.hist, histogram file:{0}/{1}.histogram\n".format(args.outdir,args.prefix))
    flog.write("use flash combine the R1 and R2 finish, use time : {0:.5f} s\n\n\n".format(time.time()-be))
    flog.close()

def cutadapt(args,commands):
    flog = open(args.logfile,'a')
    flog.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " use cutadapter find the linker start ...\n")
    flog.flush()
    be = time.time()
    returncode,returnresult = subprocess.getstatusoutput("cutadapt -b file:{0} -n 14 --no-indels -o {1}/{2}.noLinker --info-file {1}/{2}.Linker_info --discard -O {3} {4} > {1}/{2}.stat".format(args.linkerfile,args.outdir,args.prefix,args.minlen,args.input))
    if returncode != 0:
        sys.stderr.write("[ERROR]: failed to run cutadapt : \n{0}\n".format(returnresult))
        sys.stderr.flush()
    for i in returnresult.split("\n"):
        print(i)
    flog.write("result file: Linker information file:{0}/{1}.Linker_info, remove linker fastq file:{0}/{1}.noLinker\n".format(args.outdir,args.prefix))
    flog.write("use cutadapter find the linker finish, use time : {0:.5f} s\n\n\n".format(time.time()-be))
    flog.close()

def multiLinker(args,commands):
    from ChRDbox import mergeToOneLine
    flog = open(args.logfile,'a')
    flog.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " convert to multi-line to one line start ...\n")
    flog.flush()
    be = time.time()
    
    infile = args.input
    outfile = "{0}/{1}.cut.info".format(args.outdir,args.prefix)
    if infile.endswith(".gz"):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,"r")
    with open(outfile,'w') as fout:
        flag = ""
        static = {}
        readarr = []
        nn = 0
        for line in fin:
            nn += 1
            if nn%1000000 == 0:
                flog.write("load {0} M records...\n".format( nn//1000000 ))
                flog.flush()
            temp = line.strip().split()
            if flag != "" and flag != temp[0]:
                linkercount,linkerstat = mergeToOneLine(readarr,fout)
                if linkerstat:
                    if linkercount not in static:
                        static[linkercount] = 0
                    static[linkercount] += 1
                else:
                    if 'nolinker' not in static:
                        static['nolinker'] = 0
                    static['nolinker'] += 1
                readarr = []
            flag = temp[0]
            readarr.append(temp)
        flog.write("total have {0} records...\n".format(nn))
        flog.flush()
        linkercount,linkerstat = mergeToOneLine(readarr,fout)
        if linkerstat:
            if linkercount not in static:
                static[linkercount] = 0
            static[linkercount] += 1
        else:
            if 'nolinker' not in static:
                static['nolinker'] = 0
            static['nolinker'] += 1
    for i in sorted(static.keys(), key=functools.cmp_to_key(mycmp)):
        flog.write("{0}\t{1}\n".format(i,static[i]))
    fin.close()
    flog.write("result file: combine multi-linker into one line file:{0}/{1}.cut.info\n".format(args.outdir,args.prefix))
    flog.write("convert to multi-line to one line finish, use time : {0:.5f} s\n\n\n".format(time.time()-be))
    flog.close()

def splitcombine(args,commands):
    from ChRDbox import combine2reads
    flog = open(args.logfile,'a')
    flog.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " split combine file to fastq start ...\n")
    flog.flush()
    be = time.time()
    
    infile = args.input
    RNAfile = "{0}/{1}.RNA.fq".format(args.outdir,args.prefix)
    DNAfile = "{0}/{1}.DNA.fq".format(args.outdir,args.prefix)
    if infile.endswith(".gz"):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,"r")
    static = {}

    with open(RNAfile,'w') as fRNA,open(DNAfile,'w') as fDNA:
        nn = 0
        for line in fin:
            nn += 1
            if nn%1000000 == 0:
                flog.write("load {0} M records...\n".format( nn//1000000 ))
                flog.flush()
            combine2reads(line,static,fRNA,fDNA)
    fin.close()
    flog.write("total have {0} records...\n".format(nn))
    flog.flush()
    for i in sorted(static.keys(), key=functools.cmp_to_key(mycmp)):
        flog.write("{0}\t{1}\n".format(i,static[i]))
    flog.write("marker 1_1 means: contain result\n")
    flog.write("marker 1_2 means: DNA or RNA seq length less than 18bp\n")
    flog.write("marker 1_3 means: just have one linker(no DNA/RNA seq)\n")
    flog.write("result file: RNA fastq file:{0}, DNA fastq file:{1}\n".format(RNAfile,DNAfile))
    flog.write("split combine file to fastq finish, use time : {0:.5f} s\n\n\n".format(time.time()-be))
    flog.close()

def splitnotcombine(args,commands):
    from ChRDbox import notcombine2reads
    flog = open(args.logfile,'a')
    flog.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " split notcombine file to fastq start ...\n")
    flog.flush()
    be = time.time()
    
    infile1 = args.input1
    infile2 = args.input2
    RNAfile = "{0}/{1}.RNA.fq".format(args.outdir,args.prefix)
    DNAfile = "{0}/{1}.DNA.fq".format(args.outdir,args.prefix)
    matchfile = "{0}/{1}.match.txt".format(args.outdir,args.prefix)
    if infile1.endswith(".gz"):
        fin1 = gzip.open(infile1,'rt')
    else:
        fin1 = open(infile1,"r")
    if infile2.endswith(".gz"):
        fin2 = gzip.open(infile2,'rt')
    else:
        fin2 = open(infile2,"r")
    static = {}

    with open(RNAfile,'w') as fRNA,open(DNAfile,'w') as fDNA,open(matchfile,'w') as fmatch:
        nn = 0
        for line1,line2 in zip(fin1,fin2):
            nn += 1
            if nn%1000000 == 0:
                flog.write("load {0} M records...\n".format( nn//1000000 ))
                flog.flush()
            notcombine2reads(line1,line2,static,fRNA,fDNA,fmatch)
    fin1.close()
    fin2.close()
    flog.write("total have {0} records...\n".format(nn))
    flog.flush()
    for i in sorted(static.keys(), key=functools.cmp_to_key(mycmp)):
        flog.write("{0}\t{1}\n".format(i,static[i]))
    flog.write("marker 1_0_0_0 means: R1 just have linker,R2 just have seq\n")
    flog.write("marker 1_1_0_0 means: R1 DNA or RNA length less than 18bp\n")
    flog.write("marker 1_2_0_0 means: R1 contain result\n")
    flog.write("marker 0_0_1_0 means: R2 just have linker,R1 just have seq\n")
    flog.write("marker 0_0_1_1 means: R2 DNA or RNA length less than 18bp\n")
    flog.write("marker 0_0_1_2 means: R2 contain result\n")
    flog.write("marker 1_0_1_0 means: R1 and R2 both have linker, but R1 or R2 just have linker(no seq)\n")
    flog.write("marker 1-A-1-A means: R1 and R2 both have A linker, there are two or more linker\n")
    flog.write("marker 1-A_-1-A_ means: R1 and R2 both have A_ linker, there are two or more linker\n")
    flog.write("marker 1-A-1-A_-1 means: R1 and R2 have the same linker\n")
    flog.write("marker 1-A_-1-A_ means: R1 and R2 maybe have the same linker or not\n")
    flog.write("marker 1-A_-1-A-1 means: R1 and R2 have the same linker\n")
    flog.write("marker 1-A_-1-A means: R1 and R2 maybe have the same linker or not\n")
    flog.write("result file: RNA fastq file:{0}, DNA fastq file:{1}\n".format(RNAfile,DNAfile))
    flog.write("split notcombine file to fastq finish, use time : {0:.5f} s\n\n\n".format(time.time()-be))
    flog.close()

def combine(args,commands):
    flog = open(args.logfile,'a')
    flog.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " combine multi file start ...\n")
    flog.flush()
    be = time.time()

    outfile = "{0}/{1}.fastq".format(args.outdir,args.prefix)
    fout = open(outfile,'w')
    total = 0
    for i in args.input:
        index = 0
        with open(i,'r') as fin:
            for line in fin:
                index += 1
                fout.write(line)
            flog.write("{0} load finished, total {1} records ...".format(i,index//4))
    fout.close()
    flog.write("All file combine finished.{0} records".format(total//4))
    flog.write("result file: combine file:{0}/{1}.fastq\n".format(args.outdir,args.prefix))
    flog.write("use cutadapter find the linker finish, use time : {0:.5f} s\n\n\n".format(time.time()-be))
    flog.close()

def DNAmap(args,commands):
    flog = open(args.logfile,'a')
    flog.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " combine multi file start ...\n")
    flog.flush()
    be = time.time()

    longfile = "{0}/{1}.long.fastq".format(args.outdir,args.prefix)
    shortfile = "{0}/{1}.short.fastq".format(args.outdir,args.prefix)
    with open(args.input,'r') as fin,open(longfile,'w') as flong,open(shortfile,'w') as fshort:
        index = 0
        mystr = ""
        length = 0
        countlong = 0
        countshort = 0
        for index,line in enumerate(fin):
            mystr += line
            if index%4 == 1:
                length = len(line.strip())
            elif index%4 == 3:
                if length >= args.length:
                    flong.write(mystr)
                    countlong += 1
                else:
                    fshort.write(mystr)
                    countshort += 1
                mystr = ""
    flog.write("The long DNA file contain {0} records".format(countlong))
    flog.write("The short DNA file contain {0} records".format(countshort))

    returncode,returnresult = subprocess.getstatusoutput("bwa aln -t {0} -f {1}/{2}.short.sai {3} {1}/{2}.short.fastq && bwa samse -f {1}/{2}.short.sam {3} {1}/{2}.short.sai {1}/{2}.short.fastq && bwa mem -t {0} {3} {1}/{2}.long.fastq > {1}/{2}.long.sam".format(args.threads,args.outdir,args.prefix,args.genome))
    if returncode != 0:
        sys.stderr.write("[ERROR]: failed to run flash : \n{0}\n".format(returnresult))
        sys.stderr.flush()
    
    flog.write("result file: combine file:{0}/{1}.fastq\n".format(args.outdir,args.prefix))
    flog.write("use cutadapter find the linker finish, use time : {0:.5f} s\n\n\n".format(time.time()-be))
    flog.close()

def getargs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Program: ChRD-tools
Version: 1.0.0
Usage: ChRD <command> [options]
Commands:
    flash        use flash to combine R1 and R2.
    cutadapt   convent format(SPRITE)
    multiLinker         find the colocal region
    cluster         find the cluster/Loop and cluster/Loop anchor
    visual          visualization
==========================================================================


""")
    parser.add_argument('--version', action='version',version='%(prog)s {0} by xyhuang'.format('1.0.0'))
    subparsers = parser.add_subparsers(title='subcommands',description='valid subcommands',help='')

    ## add subcommand 
    parser_runMulti = subparsers.add_parser('flash',help='combine the R1 and R2')
    parser_runMulti.add_argument("-r1", "--inputR1",required=True,help="the R1 data")
    parser_runMulti.add_argument("-r2", "--inputR2",required=True,help="the R2 data")
    parser_runMulti.add_argument("-v", "--version", action = "version", version = "%(prog)s 1.0.0 by xyhuang")
    parser_runMulti.add_argument("-t", "--threads", default = 1, help = "the threads count [default:1]")
    parser_runMulti.add_argument("-max", "--maxlength",default = 145, help = "the max length when R1 and R2 overlap")
    parser_runMulti.add_argument("-p", "--prefix", default = "output", help = "the output file prefix [default:output] ")
    parser_runMulti.add_argument("-o", "--outdir", default = "./", help = "the output dir [default:./] ")
    parser_runMulti.add_argument("-l", "--logfile", default = "flash.log", help = "the script run log file [default:flash.log] ")
    parser_runMulti.set_defaults(func=flash)

    parser_format = subparsers.add_parser('cutadapt',help='linker finding')
    parser_format.add_argument('-i', "--input", required = True, help="the input file")
    parser_format.add_argument("-o", "--outdir", default = "./", help = "the output dir [default:./]")
    parser_format.add_argument('-f', "--linkerfile", required = True ,help="the linker file")
    parser_format.add_argument("-p", "--prefix", default = "output", help = "the output file prefix [default:output]")
    parser_format.add_argument('-min', "--minlen", default=18 ,type = int ,help="the min of reads length [default:18]")
    parser_format.add_argument("-l", "--logfile", default = "cutadapt.log", help = "the script run log file [default:cutadapt.log]")
    parser_format.set_defaults(func=cutadapt)

    parser_cluster = subparsers.add_parser('multiLinker',help='Some reads have multi Linker. In cutadapt result file, they split into multi-line result')
    parser_cluster.add_argument('-i', "--input", required = True, help="the input file")
    parser_cluster.add_argument("-o", "--outdir", default = "./", help = "the output dir [default:./]")
    parser_cluster.add_argument("-p", "--prefix", default = "output", help = "the output file prefix [default:output]")
    parser_cluster.add_argument("-l", "--logfile", default = "multiLinker.log", help = "the script run log file [default:multiLinker.log]")
    parser_cluster.set_defaults(func=multiLinker)

    parser_colocal = subparsers.add_parser('splitcombine',help='split combine data to fastq')
    parser_colocal.add_argument('-i', "--input", required = True, help="the input file")
    parser_colocal.add_argument("-o", "--outdir", default = "./", help = "the output dir [default:./]")
    parser_colocal.add_argument("-p", "--prefix", default = "output", help = "the output file prefix [default:output]")
    parser_colocal.add_argument("-l", "--logfile", default = "splitcombine.log", help = "the script run log file [default:splitcombine.log]")
    parser_colocal.set_defaults(func=splitcombine)

    parser_visual = subparsers.add_parser('splitnotcombine',help='split not combine data to fastq')
    parser_visual.add_argument('-i1', "--input1", required = True, help="the input1 file")
    parser_visual.add_argument('-i2', "--input2", required = True, help="the input2 file")
    parser_visual.add_argument("-o", "--outdir", default = "./", help = "the output dir [default:./]")
    parser_visual.add_argument("-p", "--prefix", default = "output", help = "the output file prefix [default:output]")
    parser_visual.add_argument("-l", "--logfile", default = "splitnotcombine.log", help = "the script run log file [default:splitnotcombine.log]")
    parser_visual.set_defaults(func=splitnotcombine)

    parser_combine = subparsers.add_parser('combine',help='combine multi file to one file')
    parser_combine.add_argument("-i", "--input",nargs='+', required=True,help="the input file")
    parser_combine.add_argument("-o", "--outdir", default = "./", help = "the output dir [default:./]")
    parser_combine.add_argument("-p", "--prefix", default = "output", help = "the output file prefix [default:output]")
    parser_combine.add_argument("-l", "--logfile", default = "combine.log", help = "the script run log file [default:combine.log]")
    parser_combine.set_defaults(func=combine)

    parser_DNAmap = subparsers.add_parser('DNAmap',help='map the DNA seq to reference genome')
    parser_DNAmap.add_argument('-i', "--input", required = True, help="the DNA seq input")
    parser_DNAmap.add_argument("-g", "--genome", required = True , help = "the reference genome path with bwa index")
    parser_DNAmap.add_argument("-o", "--outdir", default = "./", help = "the output dir [default:./]")
    parser_DNAmap.add_argument("-p", "--prefix", default = "output", help = "the output file prefix [default:output]")
    parser_DNAmap.add_argument("-t", "--threads", default = 1, type = int, help = "the threads count [default:1]")
    parser_DNAmap.add_argument("-len", "--length", default = 70, type = int, help = "less the length:use bwa aln;great the length:use bwa mem[default:70]")
    parser_DNAmap.add_argument("-l", "--logfile", default = "combine.log", help = "the script run log file [default:combine.log]")
    parser_DNAmap.set_defaults(func=DNAmap)

    commands = sys.argv[1:]
    if ((not commands) or ((commands[0] in ['runMulti', 'SPRITEconvert', 'colocal','cluster','visual']) and len(commands) == 1)):
        commands.append('-h')
    args = parser.parse_args(commands)
    return args, commands

if __name__ == "__main__":
    args, commands = getargs()
    args.func(args, commands)