import sys
import numpy as np
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot(indat,bins_arr,chr_arr,mmin=2,mmax=98):
    pmin = np.nanpercentile(indat,mmin)
    pmin = 0
    pmax = np.nanpercentile(indat,mmax)
    print(pmin,pmax)
    plt.imshow(indat, vmin=pmin, vmax=pmax, origin = 'lower', interpolation='none', cmap=plt.cm.Reds)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xticks(rotation=90)
    plt.vlines(bins_arr[:-1], 0,bins_arr[-1], colors = "grey", linestyles = "solid",lw=0.1)
    print(bins_arr[:-1])
    plt.hlines(bins_arr[:-1], 0,bins_arr[-1], colors = "grey", linestyles = "solid",lw=0.1)
    start = bins_arr.copy()
    start.insert(0,1)
    start.pop()
    plt.yticks([int((bins_arr[i]+start[i])/2) for i in range(len(bins_arr))],[chr_arr[i] for i in range(len(chr_arr))])
    plt.xticks([int((bins_arr[i]+start[i])/2) for i in range(len(bins_arr))],[chr_arr[i] for i in range(len(chr_arr))])
    plt.colorbar()
    return plt

def main(file):
    bins_arr = [] # chr position list
    chr_arr = [] # chr name list
    flag = None
    temp = []
    
    with open(file[0],'r') as fl: # get position and chrom name
        for line in fl:
            temp = line.strip().split()
            if flag != None and flag != temp[0]:
                chr_arr.append(flag)
                bins_arr.append(int(temp[3])-1)
            flag = temp[0]
    chr_arr.append(temp[0])
    bins_arr.append(int(temp[3])-1)

    aa = np.loadtxt(file[1]).T  # transform coefficient matrix
    if min(aa[0].min(), aa[1].min()) == 0:
        N = int(max(aa[0].max(),aa[1].max()) + 1)
        counts=sparse.coo_matrix((aa[2], (aa[0],aa[1])),shape=(N,N),dtype=float)
    else:
        N = int(max(aa[0].max(),aa[1].max()))
        counts = sparse.coo_matrix((aa[2],(aa[0]-1,aa[1]-1)),shape=(N,N),dtype=float)
    
    mydat = np.array(counts.todense())

    pp = PdfPages(file[2])  # print picture on pdf file
    plotraw=plot(mydat,bins_arr,chr_arr)
    # np.savetxt("raw.txt",mydat)
    plotraw.savefig(pp, format='pdf')
    plotraw.close()

    nr = mydat/(np.sqrt(mydat.sum(axis=0))*np.sqrt(mydat.sum(axis=1).reshape(-1,1)))  #  VC_sqrt normalization
    where_are_nan = np.isnan(nr)
    nr[where_are_nan] = 0
    cg = (mydat.sum()/2)/(nr.sum()/2)
    normal = nr*cg
    
    plotnor=plot(normal,bins_arr,chr_arr,30)
    plotnor.savefig(pp, format='pdf')
    plotnor.close()
    np.savetxt("{0}.nor.txt".format(file[1]),normal)
    pp.close()

if __name__ == "__main__":
    main(sys.argv[1:])
