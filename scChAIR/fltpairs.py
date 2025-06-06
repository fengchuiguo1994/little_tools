import sys

import gzip
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

fin = readFile(sys.argv[1])
fout = writeFile(sys.argv[2])
for line in fin:
    if "chrM" in line or "chrY" in line:
        continue
    if line.startswith("#"):
        fout.write(line)
    else:
        tmp = line.strip().split()
        if "-M" in tmp[1] or "-P" in tmp[1]:
            chrom1, flag1 = tmp[1].split("-")
        else:
            chrom1 = tmp[1]
            flag1 = "None"
        if "-M" in tmp[3] or "-P" in tmp[3]:
            chrom2, flag2 = tmp[3].split("-")
        else:
            chrom2 = tmp[3]
            flag2 = "None"
        if flag1 == "None" and flag2 == "None":
            continue
        elif flag1 != "None" and flag2 != "None":
            if chrom1 == chrom2 and flag1 != flag2 and abs(int(tmp[2]) - int(tmp[4])) < 5000000:
                continue
            fout.write(line)
        elif flag1 != "None":
            tmp[3] = "{0}-{1}".format(tmp[3],flag1)
            fout.write("{0}\n".format("\t".join(tmp)))
        elif flag2 != "None":
            tmp[1] = "{0}-{1}".format(tmp[1],flag2)
            fout.write("{0}\n".format("\t".join(tmp)))
        else:
            print("error")
            sys.exit(1)
fin.close()
fout.close()