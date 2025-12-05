import cooler, argparse, matplotlib
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from collections import defaultdict

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

parser = argparse.ArgumentParser(description='Plot interaction matrix from .cool file')
parser.add_argument('--cool', type=str, required=True, help='Path to the .cool file')
parser.add_argument('--output' ,'-o', type=str, required=True, help='Output figure file path')
parser.add_argument('--all', '-a', action='store_true', help='Plot whole genome interaction matrix')
parser.add_argument('--chrom', '-c', type=str, help='Chromosome name to plot if --all is not set')
parser.add_argument('--start', '-s', type=int, default=0, help='Start position (bp) if --all is not set')
parser.add_argument('--end', '-e', type=int, default=0, help='End position (bp) if --all is not set')
parser.add_argument('--balance', '-b', type=str, default='raw', help='Method to balance the matrix (raw, ICE, VC_SQRT...), corresponding column should be present in the .cool file')
parser.add_argument('--chrom2', '-c2', type=str, help='Second chromosome name to plot (for inter-chromosomal matrix), if not set, intra-chromosomal matrix will be plotted')
parser.add_argument('--start2', '-s2', type=int, default=0, help='Start position (bp) of second chromosome if --chrom2 is set')
parser.add_argument('--end2', '-e2', type=int, default=0, help='End position (bp) of second chromosome if --chrom2 is set')

args = parser.parse_args()

raw_cmap = LinearSegmentedColormap.from_list('interaction',
                ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])

def print_coordinate(pos):
    
    if pos % 1000000 == 0:
        print(pos)
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
boundary = defaultdict(tuple)
# for c in [str(i) for i in range(1,20)] + ['X','Y']:
for c in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
    boundary[c] = clr1.extent(c)

wc_label = False
if not args.all:
    chromosome = args.chrom
    left = int(args.start)
    right = int(args.end)
    if not chromosome:
        raise ValueError("When --all is not set, --chrom must be provided.")
    if not (left and right):
        print('The plot will use the whole chromosome: {}'.format(chromosome))
        region = (chromosome,)
    else:
        region = (chromosome, left, right)
        print('start at chr1:{}'.format(left))
        print('end at chr1:{}'.format(right))
    if args.chrom2:
        chromosome2 = args.chrom2
        left2 = int(args.start2)
        right2 = int(args.end2)
        if not (left2 and right2):
            print('The plot will use the whole chromosome: {}'.format(chromosome2))
            region2 = (chromosome2,)
        else:
            region2 = (chromosome2, left2, right2)
            print('start at chr2:{}'.format(left2))
            print('end at chr2:{}'.format(right2))
        if len(region) == 1 and len(region2) == 1:
            print('Alos plotting inter-chromosomal matrix for {} and {}'.format(chromosome, chromosome2))
            Mchr1 = clr1.matrix(balance=balance).fetch(region[0])
            Mchr2 = clr1.matrix(balance=balance).fetch(region2[0])
            Mcross = clr1.matrix(balance=balance).fetch(region[0], region2[0])
            Mzero = np.zeros(Mcross.shape)
            Mhead = np.hstack((Mchr1, Mzero))
            Mtail = np.hstack((Mcross.T, Mchr2))
            M1 = np.vstack((Mhead, Mtail))
            wc_label=True
        else:
            M1 = clr1.matrix(balance=balance).fetch(region, region2)
    else:
        if len(region) == 1:
            M1 = clr1.matrix(balance=balance).fetch(region[0])
        else:
            M1 = clr1.matrix(balance=balance).fetch(region)


else:
    Mtmp = clr1.matrix(balance=balance)
    M1 = Mtmp[:]

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

if wc_label:
    fontsize=6.5
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    offset = 0.02 * (xmax - xmin)
    # square1 = patches.Rectangle((xmin + 2.8*offset, ymin - 2.5*offset), 2*offset, 2*offset, edgecolor = 'none', facecolor = 'red')
    # ax.add_patch(square1)
    # ax.text(xmin + 5.3*offset, ymin - 1.5*offset, print_max(vmax1), ha = 'left', va = 'center', fontsize = fontsize)
    square1 = patches.Rectangle((xmax - 3.8*offset, ymax + 2.5*offset), 2*offset, 2*offset, edgecolor = 'none', facecolor = 'red')
    ax.add_patch(square1)
    ax.text(xmax - 5.3*offset, ymax + 3.5*offset, print_max(vmax1), ha = 'right', va = 'center', fontsize = fontsize)
    idx = (boundary[chromosome][1]*2 - 1) / 2
    ax.axvline(x=idx, color='black', linewidth=0.5)
    ax.axhline(y=idx, color='black', linewidth=0.5)
    for spine in ['right', 'top', 'bottom', 'left']:
        ax.spines[spine].set_visible(False)
    ax.hlines(y=0, xmin=0, xmax=idx, color='black', linewidth=0.5, antialiased=False, zorder=3)
    ax.vlines(x=xmax-0.5, ymin=ymin, ymax=idx-0.5, color='black', linewidth=0.5, antialiased=False, zorder=3)
    ax.hlines(y=ymin-0.5, xmin=0, xmax=xmax-0.5, color='black', linewidth=0.5, antialiased=False, zorder=3)
    ax.vlines(x=0, ymin=ymin, ymax=ymax + 0.5, color='black', linewidth=0.5, antialiased=False, zorder=3)
    label_pos = [(boundary[chromosome][1] - boundary[chromosome][0])/2, (boundary[chromosome2][1] - boundary[chromosome2][0])/2 + boundary[chromosome][1]]
    ax.set_xticks(label_pos, ['chr' + chromosome, 'chr' + chromosome2], fontsize=fontsize)
    ax.set_yticks(label_pos, ['chr' + chromosome, 'chr' + chromosome2], fontsize=fontsize)
    ax.tick_params(axis='both', which='both', length=0)


elif not args.all:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    fontsize=6.5
    offset = 0.02 * (xmax - xmin)
    ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    if len(region) != 1:
        ax.text(xmin, ymin+0.8*offset, print_coordinate(region[1]), va='top', ha='left', fontsize=fontsize)
        ax.text(xmax, ymin+0.8*offset, print_coordinate(region[2]), va='top', ha='right', fontsize=fontsize)
        ax.text(xmax+3.5*offset, ymax, print_coordinate(region[1]), rotation=270, va='top', ha='right', fontsize=fontsize)
        ax.text(xmax+3.5*offset, ymin, print_coordinate(region[2]), rotation=270, va='bottom', ha='right', fontsize=fontsize)
    else:
        pass
    ax.text((xmin + xmax)/2, ymin + 0.8*offset, region[0], va='top', ha='center', fontsize=fontsize)
    ax.text(xmax + 3.5*offset, (ymin + ymax)/2, region[0], rotation=270, va='center', ha='right', fontsize=fontsize)
    square1 = patches.Rectangle((xmax - 3.8*offset, ymax + 2.5*offset), 2*offset, 2*offset, edgecolor = 'none', facecolor = 'red')
    ax.add_patch(square1)
    ax.text(xmax - 5.3*offset, ymax + 3.5*offset, print_max(vmax1), ha = 'right', va = 'center', fontsize = fontsize)

else:
    fontsize=3
    label_pos = []
    # for i in [str(i) for i in range(1,20)] + ['X']:
    for i in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
        idx = (boundary[i][1]*2 - 1) / 2
        ax.axvline(x=idx, color='black', linewidth=0.2)
        ax.axhline(y=idx, color='black', linewidth=0.2)
        label_pos.append((boundary[i][0] + boundary[i][1]) / 2)
    # label_pos.append((boundary['Y'][0] + boundary['Y'][1]) / 2)
    # ax.set_xticks(label_pos, ['chr1'] + [str(i) for i in range(2,20)] + ['X','Y'], fontsize=fontsize)
    # ax.set_yticks(label_pos, ['chr1'] + [str(i) for i in range(2,20)] + ['X','Y'], fontsize=fontsize)
    ax.set_xticks(label_pos, ['chr01'] + ["02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"], fontsize=fontsize)
    ax.set_yticks(label_pos, ['chr01'] + ["02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"], fontsize=fontsize)
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='both', which='both', length=0)
    cax = fig.add_axes([ax.get_position().x1 + 0.03, ax.get_position().y0, 0.02, ax.get_position().height])
    cbar = fig.colorbar(sc, cax=cax, fraction=1.0)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(['0', print_max(vmax1)])
    cbar.ax.tick_params(labelsize=3, width=0.2, length=1)
    cbar.outline.set_linewidth(0.2)

plt.savefig(outfig, dpi=1600, bbox_inches='tight')
plt.close()

