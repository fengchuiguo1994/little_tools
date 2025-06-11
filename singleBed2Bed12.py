import sys

def tobed12(mylist,out):
    out.write("{0}\t{1}\t{2}\t{3}\t0\t+\t{4}\t{5}\t0\t1\t{6}\t0,\n".format(mylist[0][0],mylist[0][1],mylist[0][2],mylist[0][3],mylist[0][1],mylist[0][2],(int(mylist[0][2])-int(mylist[0][1]))))
    
def main(bed4file,bed12file):
    with open(bed4file,'r') as bed3,open(bed12file,'w') as bed12:
        for line in bed3:
            tmp = line.strip().split()
            tobed12(tmp,bed12)

if __name__ == "__main__":
    # bed3 format input file, bed12 format output file
    main(sys.argv[1],sys.argv[2])