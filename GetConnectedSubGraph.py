import networkx as nx
import sys

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

fin = readFile(sys.argv[1])

edge = []
for line in fin:
    tmp = line.strip().split()
    edge.append((tmp[0],tmp[1]))
fin.close()

G = nx.Graph()
G.add_edges_from(edge)
fout = writeFile(sys.argv[2])
for c in nx.connected_components(G):
    fout.write("{0}\n".format("\t".join(sorted(c))))
fout.close()
