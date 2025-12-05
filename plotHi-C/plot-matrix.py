import sys, cooler
import numpy as np 
import matplotlib.pylab as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import os, subprocess
import cooltools.lib.plotting
from packaging import version
import cooltools
import bioframe, pyBigWig
from IPython import display
from matplotlib.transforms import Affine2D
import matplotlib as mpl 
mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['ps.fonttype'] = 42

raw_cmap = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])

def print_coordinate(pos):
    if pos % 1000000 == 0:
        return '{0}M'.format(pos//1000000)
    else:
        return '{0:.2f}M'.format(pos/1000000)

def print_max(val):
    return '{0:.1f}'.format(val)

def print_num(num):
    num = float(num)
    n = 0
    while num > 0:
        num = num//10
        n+=1
    return n

cool_uri = sys.argv[1] # "/data/home/ruanlab/huangxingyu/MH63RIC/tissueRICresult/figure2new/2.tohic/paniclenorRNA/paniclenorRNA.genome.chromosome.4dn.hic.mcool"
uri = sys.argv[2] # 1000
outfig = sys.argv[3] # "pdf.pdf"
chromosome = sys.argv[4] # "Chr05" # "05"
left = int(sys.argv[5]) # 22921000
right =int(sys.argv[6]) # 23060000

print('start at:{}'.format(left))
print('end at:{}'.format(right))
region = (chromosome, left, right)

clr1 = cooler.Cooler(f'{cool_uri}::/resolutions/{uri}')
# M1 = clr1.matrix(balance='VC_SQRT').fetch(region)
M1 = clr1.matrix(balance=False).fetch(region)
# vmax1 = np.percentile(M1[M1>0], 90)
vmax1 = 3
print('M1 - maximum value in colorbar: {0}'.format(vmax1))

M1 = M1 / vmax1
M1[M1>1] = 1
M = M1

label_size = 6
label_pad = 1

fig = plt.figure(figsize=(1.6,1.6))
ax = fig.add_subplot()
sc = ax.imshow(M, interpolation='none', cmap=raw_cmap, aspect = 'auto', vmax=1)

ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

for spine in ['right', 'top', 'bottom', 'left']:
    ax.spines[spine].set_linewidth(0.5)

xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
fontsize=6.5
offset = 0.02 * (xmax - xmin)
ax.text(xmin, ymin+0.8*offset, print_coordinate(region[1]), va='top', ha='left', fontsize=fontsize)
ax.text(xmax, ymin+0.8*offset, print_coordinate(region[2]), va='top', ha='right', fontsize=fontsize)
ax.text(xmax+3.5*offset, ymax, print_coordinate(region[1]), rotation=270, va='top', ha='right', fontsize=fontsize)
ax.text(xmax+3.5*offset, ymin, print_coordinate(region[2]), rotation=270, va='bottom', ha='right', fontsize=fontsize)
ax.text((xmin + xmax)/2, ymin + 0.8*offset, region[0], va='top', ha='center', fontsize=fontsize)
ax.text(xmax + 3.5*offset, (ymin + ymax)/2, region[0], rotation=270, va='center', ha='right', fontsize=fontsize)
import matplotlib.patches as patches
square1 = patches.Rectangle((xmax - 3.8*offset, ymax + 2.5*offset), 2*offset, 2*offset, edgecolor = 'none', facecolor = 'red')
ax.add_patch(square1)

ax.text(xmax - 5.3*offset, ymax + 3.5*offset, print_max(vmax1), ha = 'right', va = 'center', fontsize = fontsize)

plt.savefig(outfig, dpi=1600, bbox_inches='tight')
plt.close()
