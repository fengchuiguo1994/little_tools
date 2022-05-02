import sys
import re

def main(file):
    hpv = {}
    with open(file,'r') as fl:
        for line in fl:
            line = line.strip()
            line2 = fl.next()
            genome = (re.findall(r"<title>(.+?),",line2))[0]
            hpv.setdefault(genome,[]).append(line)
    for i,j in hpv.items():
        print(i,j)

if __name__ == "__main__":
    main(sys.argv[1])