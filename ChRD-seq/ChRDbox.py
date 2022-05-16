import pysam
import time
import re
import os
import subprocess
import random
from operator import itemgetter
from bx.intervals.intersection import Intersecter, Interval
import sys
import pybedtools
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import matplotlib.path as mpath
import pandas as pd
import functools
'''
function : give runMulti the deal function
create : huang xingyu
version : 1.0.1
begin : 2019-07-11
end : 
'''
def mergeToOneLine(mylist,fout):
    num = len(mylist)
    if num == 1: # just have one linker or no linker
        temp = mylist[0]
        if len(temp) == 5: # no linker
            fout.write("{0}\t0\t{1}\t{2}\n".format(temp[0],temp[3],temp[4]))
            return 1,False
        elif len(temp) == 12: #------ -- ------
            fout.write("{0}\t1\t{1}\t{2}#{3}\t{4}\t{5}\t{6}\t{7}\n".format(temp[0],temp[5],temp[6],temp[8],temp[7],temp[9],temp[10],temp[11]))
        elif len(temp) == 10 and temp[3] == "0": # -- ------
            fout.write("{0}\t1\t#\t{1}#{2}\t{3}\t#\t{4}\t{5}\n".format(temp[0],temp[5],temp[7],temp[6],temp[8],temp[9]))
        elif len(temp) == 10 and temp[3] != "0": # ------ --
            fout.write("{0}\t1\t{1}\t{2}#{3}\t#\t{4}\t{5}\t#\n".format(temp[0],temp[5],temp[6],temp[7],temp[8],temp[9]))
        elif len(temp) == 8: # -- (only have linker)
            fout.write("{0}\t1\t{1}#{2}\t{3}\n".format(temp[0],temp[5],temp[6],temp[7]))
        elif len(temp) == 3 or len(temp) == 2: # no seq or no linker
            fout.write("{0}\t0\t#\t#\n".format(temp[0]))
        else:
            raise Exception("it maybe have some error1",temp)
        return 1,True
    elif num > 1:
        out = []
        value = []
        for i in range(num):
            temp = mylist[i]
            if i == 0:
                if len(temp) == 12: # ------ -- ------
                    temp[6] += "#{0}".format(temp[8])
                    out.extend([temp[5],temp[6],temp[7]])
                    value.extend([temp[9],temp[10],temp[11]])
                elif len(temp) == 10 and temp[3] == "0": # -- ------
                    temp[5] += "#{0}".format(temp[7])
                    out.extend([temp[5],temp[6]])
                    value.extend([temp[8],temp[9]])
                elif len(temp) == 10 and temp[3] != "0": # ------ --
                    temp[6] += "#{0}".format(temp[7])
                    out.extend([temp[5],temp[6]])
                    value.extend([temp[8],temp[9]])
                elif len(temp) == 8: # -- (just have linker)
                    temp[5] += "#{0}".format(temp[6])
                    out.extend(temp[5])
                    value.extend(temp[7])
                else:
                    raise Exception("it maybe have some error2",temp)
            else:
                if len(temp) == 12: # ------ -- ------
                    str = "{0}{1}{2}".format(temp[5],temp[6],temp[7])
                    temp[6] += "#{0}".format(temp[8])
                    for j in range(len(out)):
                        if out[j] == str:
                            arr = out[j+1:]
                            out = out[0:j]
                            out.extend([temp[5],temp[6],temp[7]])
                            out.extend(arr)
                            arr = value[j+1:]
                            value = value[0:j]
                            value.extend([temp[9],temp[10],temp[11]])
                            value.extend(arr)
                            break
                elif len(temp) == 10 and temp[3] == "0": # -- ------
                    str = "{0}{1}".format(temp[5],temp[6])
                    temp[5] += "#{0}".format(temp[7])
                    for j in range(len(out)):
                        if out[j] == str:
                            arr = out[j+1:]
                            out = out[0:j]
                            out.extend([temp[5],temp[6]])
                            out.extend(arr)
                            arr = value[j+1:]
                            value = value[0:j]
                            value.extend([temp[8],temp[9]])
                            value.extend(arr)
                            break
                elif len(temp) == 10 and temp[3] != "0": # ------ --
                    str = "{0}{1}".format(temp[5],temp[6])
                    temp[6] += "#{0}".format(temp[7])
                    for j in range(len(out)):
                        if out[j] == str:
                            arr = out[j+1:]
                            out = out[0:j]
                            out.extend([temp[5],temp[6]])
                            out.extend(arr)
                            arr = value[j+1:]
                            value = value[0:j]
                            value.extend([temp[8],temp[9]])
                            value.extend(arr)
                            break
                elif len(temp) == 8: # -- only linker
                    str = temp[5]
                    temp[5] += "#{0}".format(temp[6])
                    for j in range(len(out)):
                        if out[j] == str:
                            out[j] = temp[5]
                            value[j] = temp[7]
                            break
                else:
                    raise Exception("it maybe have some error2",temp)
        fout.write("{0}\t{1}\t{2}\t{3}\n".format(mylist[0][0],num,"\t".join(out),"\t".join(value)))
        return num,True
    else:
        raise Exception("it maybe have some error4",mylist)

def rc(seq):
    aa = {'A':'T','T':'A','C':'G','G':'C','a':'T','t':'A','c':'G','g':'C'}
    result = ""
    le = len(seq)
    for i in range(le):
        if seq[le-1-i] in aa:
            result += aa[seq[le-1-i]]
        else:
            result += seq[le-1-i]
    return result

def printReads(mark,temp,RNAf,DNAf):
    if mark == "A": # RNA kmer DNA
        temp[4] = rc(temp[4])
        temp[7] = temp[7][::-1]
        RNAf.write("@{0} 1\n".format(temp[0]))
        RNAf.write("{0}\n".format(temp[2]))
        RNAf.write("+\n")
        RNAf.write("{0}\n".format(temp[5]))
        DNAf.write("@{0} 2\n".format(temp[0]))
        DNAf.write("{0}\n".format(temp[4]))
        DNAf.write("+\n")
        DNAf.write("{0}\n".format(temp[7]))
    elif mark == "A_": # DNA kmer RNA
        temp[4] = rc(temp[4])
        temp[7] = temp[7][::-1]
        RNAf.write("@{0} 1\n".format(temp[0]))
        RNAf.write("{0}\n".format(temp[4]))
        RNAf.write("+\n")
        RNAf.write("{0}\n".format(temp[7]))
        DNAf.write("@{0} 2\n".format(temp[0]))
        DNAf.write("{0}\n".format(temp[2]))
        DNAf.write("+\n")
        DNAf.write("{0}\n".format(temp[5]))
    else:
        raise Exception("the result maybe have bug,error1",temp)

def combine2reads(line,mydict,RNAf,DNAf):
    temp = line.strip().split()
    num = int(temp[1])
    if num not in mydict:
        mydict[num] = 0
    mydict[num] += 1
    if num == 1: # one linker
        if len(temp) == 4: # -- just have one linker(no DNA/RNA seq)
            if '1_3' not in mydict:
                mydict['1_3'] = 0
            mydict['1_3'] += 1
        elif len(temp[2]) >= 18 and len(temp[4]) >= 18:
            if re.search(r'^\w+#(A_?)$',temp[3]):
                mark = re.search(r'^\w+#(A_?)$',temp[3])[1]
                printReads(mark,temp,RNAf,DNAf)
            if '1_1' not in mydict:
                mydict['1_1'] = 0
            mydict['1_1'] += 1
        else:
            if '1_2' not in mydict:
                mydict['1_2'] = 0
            mydict['1_2'] += 1

def compseq(s1,s2):
    total = 0
    for i,j in zip(s1,s2):
        if i == j:
            total += 1
    if total/len(s1) >= 0.8:
        return 1
    else:
        return -1

def compare(list1,list2,RNAf,DNAf):
    flag1 = 0
    flag2 = 0
    flag = 0
    if len(list1[2]) >= 10 and len(list2[4]) >= 10: # -- ------ and ------ -- ------
        tmplen = len(list1[2])
        if len(list2[4]) < len(list1[2]):
            tmplen = len(list2[4])
        tmpstr1 = list1[2][len(list1[2])-tmplen:] # R1 seq near linker
        tmpstr1 = rc(tmpstr1)
        tmpstr2 = list2[4][0:tmplen] # R2 seq near linker
        flag1 = compseq(tmpstr1,tmpstr2)

    if len(list1[4]) >= 10 and len(list2[2]) >= 10:
        tmplen = len(list1[4])
        if len(list2[2]) < len(list1[4]):
            tmplen = len(list2[2])
        tmpstr1 = list2[2][len(list2[2])-tmplen:] # R2 seq near linker
        tmpstr1 = rc(tmpstr1)
        tmpstr2 = list1[4][0:tmplen] # R1
        flag2 = compseq(tmpstr1,tmpstr2)

    if (flag1 == 1 or flag2 == 1) and (flag1 != -1 and flag2 != -1):
        rna1 = "@{0} 1\n".format(list1[0])
        rna2 = "{0}\n".format(list1[2])
        rna3 = "+\n"
        rna4 = "{0}\n".format(list1[5])
        if len(list1[2]) < len(list2[4]):
            list2[4] = rc(list2[4])
            list2[7] = list2[7][::-1]
            rna1 = "@{0} 1\n".format(list2[0])
            rna2 = "{0}\n".format(list2[4])
            rna3 = "+\n"
            rna4 = "{0}\n".format(list2[7])
        list1[4] = rc(list1[4])
        list1[7] = list1[7][::-1]
        dna1 = "@{0} 2\n".format(list1[0])
        dna2 = "{0}\n".format(list1[4])
        dna3 = "+\n"
        dna4 = "{0}\n".format(list1[7])
        if len(list2[2]) > len(list1[4]):
            dna1 = "@{0} 2\n".format(list2[0])
            dna2 = "{0}\n".format(list2[2])
            dna3 = "+\n"
            dna4 = "{0}\n".format(list2[5])
        if len(dna2) > 18 and len(rna2) > 18:
            DNAf.write("{0}{1}{2}{3}".format(dna1,dna2,dna3,dna4))
            RNAf.write("{0}{1}{2}{3}".format(rna1,rna2,rna3,rna4))
            flag = 1
    return flag,flag1,flag2

def notcombine2reads(line1,line2,mydict,RNAf,DNAf,fmatch):
    temp1 = line1.strip().split()
    temp2 = line2.strip().split()
    num1 = int(temp1[1])
    num2 = int(temp2[1])
    mark = "{0}_{1}".format(num1,num2)
    if mark not in mydict:
        mydict[mark] = 0
    mydict[mark] += 1
    total = num1+num2
    if total not in mydict:
        mydict[total] = 0
    mydict[total] += 1

    if num1 == 1 and num2 == 0: # R1 one linker , R2 no linker
        if len(temp1) == 4 and len(temp2) == 4: # R1 just have linker,R2 just have seq
            if "1_0_0_0" not in mydict:
                mydict["1_0_0_0"] = 0
            mydict["1_0_0_0"] += 1
        elif len(temp1[2]) < 18 or len(temp1[4]) < 18: # DNA/RNA too short
            if "1_1_0_0" not in mydict:
                mydict["1_1_0_0"] = 0
            mydict["1_1_0_0"] += 1
        elif re.search(r'^\w+#(A_?)$',temp1[3]): # ------ -- ------ ok
            mark = re.search(r'^\w+#(A_?)$',temp1[3])[1]
            printReads(mark,temp1,RNAf,DNAf)
            if "1_2_0_0" not in mydict:
                mydict["1_2_0_0"] = 0
            mydict["1_2_0_0"] += 1
        else:
            raise Exception("the result maybe have bug,error2",temp1,temp2)
    elif num1 == 0 and num2 == 1: # R1 no linker, R2 one linker
        if len(temp1) == 4 and len(temp2) == 4: # R1 have seq,R2 have linker
            if "0_0_1_0" not in mydict:
                mydict["0_0_1_0"] = 0
            mydict["0_0_1_0"] += 1
        elif len(temp2[2]) < 18 or len(temp2[4]) < 18: # DNA/RNA too short
            if "0_0_1_1" not in mydict:
                mydict["0_0_1_1"] = 0
            mydict["0_0_1_1"] += 1
        elif re.search(r'^\w+#(A_?)$',temp2[3]): # ------ -- ------ ok
            mark = re.search(r'^\w+#(A_?)$',temp2[3])[1]
            printReads(mark,temp2,RNAf,DNAf)
            if "0_0_1_2" not in mydict:
                mydict["0_0_1_2"] = 0
            mydict["0_0_1_2"] += 1
        else:
            raise Exception("the result maybe have bug,error2",temp1,temp2)
    elif num1 == 1 and num2 == 1: # R1 and R2 both have one linker.there are two case
        if len(temp1) == 4 or len(temp2) == 4: # R1 or R2 just have linker
            if "1_0_1_0" not in mydict:
                mydict["1_0_1_0"] = 0
            mydict["1_0_1_0"] += 1
        else:
            mark1 = re.search(r'^\w+#(A_?)$',temp1[3])[1]
            mark2 = re.search(r'^\w+#(A_?)$',temp2[3])[1]
            if mark1 == "A" and mark2 == "A": # drop (have two linker or more)
                if "1-A-1-A" not in mydict:
                    mydict["1-A-1-A"] = 0
                mydict["1-A-1-A"] += 1
            elif mark1 == "A_" and mark2 == "A_": # drop(same above)
                if "1-A_-1-A_" not in mydict:
                    mydict["1-A_-1-A_"] = 0
                mydict["1-A_-1-A_"] += 1
            elif mark1 == "A" and mark2 == "A_":
                flag,flag1,flag2 = compare(temp1,temp2,RNAf,DNAf)
                fmatch.write("{0}\t{1}\n".format(flag1,flag2))
                if flag == 1:
                    if "1-A-1-A_-1" not in mydict:
                        mydict["1-A-1-A_-1"] = 0
                    mydict["1-A-1-A_-1"] += 1
                if "1-A-1-A_" not in mydict:
                    mydict["1-A-1-A_"] = 0
                mydict["1-A-1-A_"] += 1
            elif mark1 == "A_" and mark2 == "A":
                flag,flag1,flag2 = compare(temp2,temp1,RNAf,DNAf)
                fmatch.write("{0}\t{1}\n".format(flag1,flag2))
                if flag == 1:
                    if "1-A_-1-A-1" not in mydict:
                        mydict["1-A_-1-A-1"] = 0
                    mydict["1-A_-1-A-1"] += 1
                if "1-A_-1-A" not in mydict:
                    mydict["1-A_-1-A"] = 0
                mydict["1-A_-1-A"] += 1
            else:
                raise Exception("the result maybe have bug,error3",temp1,temp2)