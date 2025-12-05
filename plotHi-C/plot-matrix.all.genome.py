import cooler, argparse, cooltools, os, subprocess
import numpy as np 
import matplotlib.pylab as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import cooltools.lib.plotting
from packaging import version
import bioframe, pyBigWig
from IPython import display
from matplotlib.transforms import Affine2D
import matplotlib.patches as patches
from collections import defaultdict
import matplotlib as mpl 
mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['ps.fonttype'] = 42

parser = argparse.ArgumentParser(description='Plot interaction matrix from .cool file')
parser.add_argument('--cool', type=str, required=True, help='Path to the .cool file')
parser.add_argument('--output' ,'-o', type=str, required=True, help='Output figure file path')
parser.add_argument('--all', '-a', action='store_true', help='Plot whole genome interaction matrix')
parser.add_argument('--chrom', '-c', type=str, help='Chromosome name to plot if --all is not set')
parser.add_argument('--start', '-s', type=int, default=0, help='Start position (bp) if --all is not set')
parser.add_argument('--end', '-e', type=int, default=0, help='End position (bp) if --all is not set')
parser.add_argument('--balance', '-b', type=str, default='raw', help='Method to balance the matrix (raw, ICE, VC_SQRT...), corresponding column should be present in the .cool file')

args = parser.parse_args()

raw_cmap = LinearSegmentedColormap.from_list('interaction',
                ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])

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

cool_uri = args.cool
outfig = args.output
balance = args.balance
if balance == 'raw':
    balance = False

clr1 = cooler.Cooler(cool_uri)

if not args.all:
    chromosome = args.chrom
    left = int(args.start)
    right = int(args.end)
    if not (chromosome and left and right):
        raise ValueError("When --all is not set, --chrom, --start, and --end must be provided.")
    region = (chromosome, left, right)
    print('start at:{}'.format(left))
    print('end at:{}'.format(right))
    M1 = clr1.matrix(balance=balance).fetch(region)

else:
    Mtmp = clr1.matrix(balance=balance)
    M1 = Mtmp[:]
    boundary = defaultdict(tuple)
    for c in [str(i) for i in range(1,20)] + ['X','Y']:
        boundary[c] = clr1.extent(c)

vmax1 = np.percentile(M1[M1>0], 95)
print('M1 - maximum value in colorbar: {0}'.format(vmax1))

M1 = M1 / vmax1
M1[M1>1] = 1
M = M1

label_size = 6
label_pad = 1

fig = plt.figure(figsize=(2,2))
ax = fig.add_subplot()
sc = ax.imshow(M, interpolation='none', cmap=raw_cmap, aspect = 'auto', vmax=1)


for spine in ['right', 'top', 'bottom', 'left']:
    ax.spines[spine].set_linewidth(0.5)


if not args.all:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    fontsize=6.5
    offset = 0.02 * (xmax - xmin)
    ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    ax.text(xmin, ymin+0.8*offset, print_coordinate(region[1]), va='top', ha='left', fontsize=fontsize)
    ax.text(xmax, ymin+0.8*offset, print_coordinate(region[2]), va='top', ha='right', fontsize=fontsize)
    ax.text(xmax+3.5*offset, ymax, print_coordinate(region[1]), rotation=270, va='top', ha='right', fontsize=fontsize)
    ax.text(xmax+3.5*offset, ymin, print_coordinate(region[2]), rotation=270, va='bottom', ha='right', fontsize=fontsize)
    ax.text((xmin + xmax)/2, ymin + 0.8*offset, region[0], va='top', ha='center', fontsize=fontsize)
    ax.text(xmax + 3.5*offset, (ymin + ymax)/2, region[0], rotation=270, va='center', ha='right', fontsize=fontsize)
    square1 = patches.Rectangle((xmax - 3.8*offset, ymax + 2.5*offset), 2*offset, 2*offset, edgecolor = 'none', facecolor = 'red')
    ax.add_patch(square1)
    ax.text(xmax - 5.3*offset, ymax + 3.5*offset, print_max(vmax1), ha = 'right', va = 'center', fontsize = fontsize)

else:
    fontsize=3
    label_pos = []
    for i in [str(i) for i in range(1,20)] + ['X']:
        idx = (boundary[i][1]*2 - 1) / 2
        ax.axvline(x=idx, color='black', linewidth=0.2)
        ax.axhline(y=idx, color='black', linewidth=0.2)
        label_pos.append((boundary[i][0] + boundary[i][1]) / 2)
    label_pos.append((boundary['Y'][0] + boundary['Y'][1]) / 2)
    ax.set_xticks(label_pos, ['chr1'] + [str(i) for i in range(2,20)] + ['X','Y'], fontsize=fontsize)
    ax.set_yticks(label_pos, ['chr1'] + [str(i) for i in range(2,20)] + ['X','Y'], fontsize=fontsize)
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='both', which='both', length=0)
    cax = fig.add_axes([ax.get_position().x1 + 0.03, ax.get_position().y0, 0.02, ax.get_position().height])
    cbar = fig.colorbar(sc, cax=cax, fraction=1.0)
    cbar.ax.tick_params(labelsize=3, width=0.2, length=1)
    cbar.outline.set_linewidth(0.2)

plt.savefig(outfig, dpi=1600, bbox_inches='tight')
plt.close()
