import sys

### Operating common files  ###
import gzip
def readFile(infile):
    """
    infile: input file
    return: file handle
    """
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
        
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

def deal(flag,count,fout):
    if len(count) == 1:
        fout.write("{0}\t{1}\n".format(flag,count[0]))
    elif len(count) > 1:
        new = []
        for i in count:
            new.append(float(i))
        fout.write("{0}\t{1}\n".format(flag,sum(new)/len(new)))
    return len(count)

fin = readFile(sys.argv[1])
fout = writeFile(sys.argv[2])

flag = None
count = []
mydict = {}
for line in fin:
    tmp = line.strip().split()
    uid = "\t".join(tmp[:3])
    if flag != None and flag != uid:
        nn = deal(flag,count,fout)
        if nn not in mydict:
            mydict[nn] = 0
        mydict[nn] += 1
        count = []
    flag = uid
    count.append(tmp[3])
nn = deal(flag,count,fout)
if nn not in mydict:
    mydict[nn] = 0
mydict[nn] += 1

for k in sorted(mydict.keys()):
    print("{0}\t{1}".format(k,mydict[k]))

fin.close()
fout.close()
sys.stdout.close()