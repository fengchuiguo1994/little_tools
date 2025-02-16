### parse the string which is the nian col in gtf file  ###
from collections import OrderedDict
def parseString(mystr):
    """
    mystr: gtf attribution information string
    return: a dictionary
    """
    d = OrderedDict()
    if mystr.endswith(";"):
        mystr = mystr[:-1]
    for i in mystr.strip().split(";"):
        tmp = i.strip().split()
        if tmp == "":
            continue
        d[tmp[0]] = tmp[1].replace("\"","")
    return d

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

### Operating sam/bam files  ###
import pysam
def readSam(insamfile):
    """
    insamfile: input sam/bam file
    return: file handle
    """
    if insamfile.endswith(".bam"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam.gz"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam"):
        insam = pysam.AlignmentFile(insamfile,'r')
    else:
        raise ValueError("the input sam/bam file is not end with sam or bam!")
    return insam
        
def writeSam(outsamfile,header):
    """
    outsamfile: output sam/bam file
    header: the sam/bam file's header(chromosome information, created by insam.handle)
    return: file handle
    """
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

class SamRead:
    """
    get sam read
    """
    def __init__(self,infile):
        self.fin = readSam(infile)
    def __iter__(self):
        flag = None
        outlist = []
        for read in self.fin:
            if flag != None and flag != read.query_name:
                yield outlist
                outlist = []
            outlist.append(read)
            flag = read.query_name
        self.fin.close()
        yield outlist
    def header(self):
        return self.fin.header

import sys
class FaRead:
    def __init__(self,infile):
        if infile == None or infile == "-":
            self.fin = sys.stdin
        else:
            self.fin = readFile(infile)
    def __iter__(self):
        flag = None
        outstr = ""
        for line in self.fin:
            if flag == None:
                flag = line.strip()[1:]
            elif line.startswith(">"):
                yield flag,outstr
                flag = line.strip()[1:]
                outstr = ""
            else:
                outstr += line.strip()
        self.fin.close()
        yield flag,outstr

class FqRead:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip().split()[0]
            rid = rid[1:]
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()

class FqRead2:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip()
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()
        
def rc(sequence):
    """
    Reverse complementary sequence
    sequence: 
    return: 
    """
    seq = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = seq.translate(trantab)
    return string

### Interval tree ###
from bx.intervals.intersection import Intersecter, Interval
def regionTree(tmp,resFrag):
    """
    tmp: The length of the list is at least 3(chrom,start,end)
    resFrag: a dictionary to store the Interval tree
    """
    if tmp[0] not in resFrag:
        resFrag[tmp[0]] = Intersecter()
    resFrag[tmp[0]].add_interval(Interval(int(tmp[1]),int(tmp[2]),tmp[3:]))

def regionTreePro(tmp,resFrag):
    """
    tmp: The length of the list is at least 3(chrom,start,end)
    resFrag: a dictionary to store the Interval tree
    """
    s = int(tmp[1])
    e = int(tmp[2])
    if tmp[0] not in resFrag:
        resFrag[tmp[0]] = Intersecter()
    result = resFrag[tmp[0]].find(s,e)
    if len(result) == 0:
        resFrag[tmp[0]].add_interval(Interval(s,e,[tmp[3:]])) # 这里用的列表，也可替换为字典、集合等等
    else:
        flag = True
        for ii in result:
            if s == ii.start and e == ii.end:
                ii.value.append(tmp[3:])
                flag = False
                break
        if flag:
            resFrag[tmp[0]].add_interval(Interval(s,e,[tmp[3:]])) # 这里用的列表，也可替换为字典、集合等等
    
def regionFind(tree,start,end):
    """
    tree: the intervals tree
    start: the given region's start
    end: the given region's end
    return: a list of overlap region
    """
    return tree.find(start,end)

### reverse strand ###
def reverseStrand(mystr):
    """
    mystr: the strand( + or -)
    return: the opposite strand(- or +)
    """
    if mystr == "+":
        return '-'
    elif mystr == "-":
        return "+"
    else:
        raise ValueError("please check the input strand, must be +/-\n{0}\n".format(mystr))

### determine whether there is overlap between the two intervals ###
def overlap(c1,s1,e1,c2,s2,e2):
    if c1 == c2:
        return e1 > s2 and s1 < e2
    return False

def detOverlap(s1,e1,s2,e2):
    """
    s1,e1,s2,e2: two intervals start and end pos
    return: the True/False whether two intervals has overlap
    """
    if s1 <= s2:
        if e1 > s2:
            return True
        else:
            return False
    else:
        if e2 > s1:
            return True
        else:
            return False
def detOverlapDis(s1,e1,s2,e2):
    """
    s1,e1,s2,e2: two intervals start and end pos
    return: the True/False whether two intervals has overlap and the overlap length or distance
    """
    if s1 <= s2:
        if e1 > s2:
            if e1 < e2:
                return True,e1-s2
            else:
                return True,e2-s2
        else:
            return False,s2-e1+1
    else:
        if e2 > s1:
            if e2 < e1:
                return True,e2-s1
            else:
                return True,e1-s1
        else:
            return False,s1-e2+1
def detOverlapChrom(c1,s1,e1,c2,s2,e2):
    """
    c1,s1,e1,c2,s2,e2: two intervals chrom start and end pos
    return: the True/False whether two intervals has overlap
    """
    if c1 == c2:
        return detOverlap(s1,e1,s2,e2)
    else:
        return False
def detOverlapDisChrom(c1,s1,e1,c2,s2,e2):
    """
    c1,c2,s1,e1,s2,e2: two intervals chrom start and end pos
    return: the True/False whether two intervals has overlap and the overlap length or distance.
    if local on two different chromosome, the distance is -1
    """
    if c1 == c2:
        return detOverlapDis(s1,e1,s2,e2)
    else:
        return False,-1

### merge the overlap interva ( like bedtools merge ) ###
def mergeRegion(regions,order=True):
    """
    regions: A list of locations [[start1,end1],[start2,end2],[start3,end3]...]
    order: is the regions in order
    return: the merge regions
    """
    import copy
    if len(regions) < 2:
        return copy.deepcopy(regions)
    if order != True:
        tmpr = sorted(regions,key=lambda x:(x[0],x[1]))
    else:
        tmpr = copy.deepcopy(regions)
    merge = []
    start = tmpr[0][0]
    end = tmpr[0][1]
    for s,e in tmpr:
        if s > end:
            merge.append([start,end])
            start = s
            end = e
        else:
            if e > end:
                end = e
    merge.append([start,end])
    return merge
### extract blank area ( like get the intron region from exon regions ) ###
def getBlankRegin(regions,order=True,merge=True):
    """
    regions: A list of locations [[start1,end1],[start2,end2],[start3,end3]...]
    order: is the regions in order
    merge: the region has no overlap
    return: the black regions
    """
    if order != True:
        raise ValueError("the input regions must be sorted")
    if merge != True:
        raise ValueError("the input regions has no overlap")
    black = []
    if len(regions) > 1:
        start = regions[0][1]
        for s,e in regions[1:]:
            black.append([start,s])
            start = e
    return black


### interval class and some method ###
class region():
    """
    A region
    start: a region start pos
    end: a region end pos
    """
    def __init__(self,start,end):
        self.start = int(start)
        self.end = int(end)
        if self.start > self.end:
            raise ValueError("start is larger than end")
    def __repr__(self):
        return "class {0}: {1}\t{2}".format(self.__class__.__name__,self.start,self.end)
    def __str__(self):
        return self.__repr__()
    def __eq__(self,other):
        if not isinstance(other,region):
            return False
        return self.start == other.start and self.end == other.end
    def __neq__(self,other):
        if not isinstance(other,region):
            return True
        return not self.__eq__(other)
    def __len__(self):
        return self.end-self.start
    def length(self):
        return self.__len__

class bed3(region):
    """
    A bed3 format region
    chrom: a chromosme
    start: a region start pos
    end: a region end pos
    """
    def __init__(self,chrom,start,end):
        region.__init__(self,start,end)
        self.chrom = chrom
    def __repr__(self):
        return "class {0}: {1}\t{2}\t{3}".format(self.__class__.__name__,self.chrom,self.start,self.end)
    def __str__(self):
        return self.__repr__()
    def __eq__(self,other):
        if not isinstance(other,bed3):
            return False
        return super(bed3, self).__eq__(other) and self.chrom == other.chrom
    def __neq__(self,other):
        if not isinstance(other,bed3):
            return True
        return not self.__eq__(other)

class bed6(bed3):
    """
    A bed6 format region
    chrom: a chromosome
    start: a region start pos
    end: a region end pos
    name: region name
    score: score
    strand: + or - strand
    """
    def __init__(self,chrom,start,end,name=".",score=0,strand="."):
        bed3.__init__(self,chrom,start,end)
        self.name = name
        self.score = score
        self.strand = strand
    def __repr__(self):
        return "class {0}: {1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(self.__class__.__name__,self.chrom,self.start,self.end,self.name,self.score,self.strand)
    def __str__(self):
        return self.__repr__()
    def __eq__(self,other):
        if not isinstance(other,bed6):
            return False
        return super(bed6, self).__eq__(other) and self.name == other.name and self.score == other.score and self.strand == other.strand
    def __neq__(self,other):
        if not isinstance(other,bed6):
            return True
        return not self.__eq__(other)
    

# color
import matplotlib.pyplot as plt

color_list = [(0.1,0.1,0.1),(0.1,0.1,0.1,0.5),'#FFFF00','#FFFF00FF','#F0F','0.5','r','blue','C0','xkcd:blue','tab:blue']
for i, j in enumerate(color_list):
    plt.plot([i,i], c=j)
    plt.annotate(repr(j), (0,i+0.1))
plt.show()


## 颜色
import numpy as np
import matplotlib.pyplot as plt

def generate_color_gradient(color1, color2, steps):
    """生成颜色渐变"""
    gradient = np.linspace(color1, color2, steps)
    return gradient

# 定义起始和结束颜色 (RGB)
color_start = np.array([1, 0, 0])  # 红色
color_end = np.array([0, 0, 1])    # 蓝色

# 生成颜色色阶
steps = 10
gradient = generate_color_gradient(color_start, color_end, steps)

# 绘制色阶
plt.figure(figsize=(8, 2))
plt.imshow([gradient], aspect='auto')
plt.axis('off')  # 隐藏坐标轴
plt.show()



## 连续型变量与颜色对应
import numpy as np
import matplotlib.pyplot as plt

def generate_color_gradient(color_start, color_end, steps):
    """生成颜色渐变"""
    return np.linspace(color_start, color_end, steps)
def map_values_to_colors(values, color_gradient):
    """将值映射到颜色上"""
    norm = plt.Normalize(np.min(values), np.max(values))
    return [color_gradient[int(norm(value) * (len(color_gradient) - 1))] for value in values]

# 设置起始和结束颜色
color_start = np.array([1, 0, 0])  # 红色
color_end = np.array([0, 0, 1])    # 蓝色
steps = 100  # 色阶分段数

# 生成颜色渐变
color_gradient = generate_color_gradient(color_start, color_end, steps)

# 示例连续数据
values = np.linspace(0, 10, 100)  # 生成从0到10的100个连续值
mapped_colors = map_values_to_colors(values, color_gradient)

# 绘制结果
plt.figure(figsize=(10, 5))
for i, color in enumerate(mapped_colors):
    plt.bar(i, 1, color=color)
plt.xticks(np.arange(0, 100, 10), np.round(np.linspace(0, 10, 10), 2))  # 设置x轴刻度
plt.yticks([])
plt.title('Continuous Data Mapped to Color Gradient')
plt.show()

## 离散型变量与颜色对应
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def generate_color_gradient(color_start, color_end, steps):
    """生成颜色渐变"""
    return np.linspace(color_start, color_end, steps)
def map_values_to_colors(values, color_gradient):
    """将值映射到颜色上"""
    norm = plt.Normalize(np.min(values), np.max(values))
    return [color_gradient[int(norm(value) * (len(color_gradient) - 1))] for value in values]

# 设置起始和结束颜色
color_start = np.array([1, 0, 0])  # 红色
color_end = np.array([0, 0, 1])    # 蓝色
steps = 100  # 色阶分段数

# 生成颜色渐变
color_gradient = generate_color_gradient(color_start, color_end, steps)

# 示例离散数据
values = np.array([2, 3, 4, 5, 6, 7, 8, 9, 2, 8, 7])
# values = np.array(["tmp{0}".format(x) for x in [2, 3, 4, 5, 6, 7, 8, 9, 2, 8, 7]])
mapped_colors = map_values_to_colors(values, color_gradient)

# 绘制结果
plt.figure(figsize=(8, 2))
for i, color in enumerate(mapped_colors):
    print(i)
    print(color)
    plt.bar(i, 1, color=color)
    
plt.xticks(values)
plt.yticks([])
plt.title('Discrete Data Mapped to Color Gradient')
plt.show()



import numpy as np
import matplotlib.pyplot as plt

def generate_color_gradient(color_start, color_end, steps):
    """生成颜色渐变"""
    return np.linspace(color_start, color_end, steps)
def map_strings_to_colors(strings, color_gradient):
    """将字符串映射到颜色上"""
    unique_strings = list(set(strings))  # 获取唯一字符串
    norm = plt.Normalize(0, len(unique_strings) - 1)
    return [color_gradient[int(norm(i) * (len(color_gradient) - 1))] for i in range(len(unique_strings))]

# 设置起始和结束颜色
color_start = np.array([1, 0, 0])  # 红色
color_end = np.array([0, 0, 1])    # 蓝色
steps = 100  # 色阶分段数

# 生成颜色渐变
color_gradient = generate_color_gradient(color_start, color_end, steps)

# 示例字符串数据
strings = ['apple', 'orange', 'banana', 'grape', 'cherry', 'apple', 'banana', 'kiwi']
mapped_colors = map_strings_to_colors(strings, color_gradient)

# 绘制结果
plt.figure(figsize=(10, 5))
for i, color in enumerate(mapped_colors):
    plt.bar(i, 1, color=color)
    
plt.xticks(np.arange(len(strings)), strings)  # 设置x轴刻度为字符串
plt.yticks([])
plt.title('String Data Mapped to Color Gradient')
plt.show()



## 颜色转换。格式转化，颜色代码到rgb，rgb转颜色代码等等。
import matplotlib.colors as mcolors

def hex_to_rgb(hex_color):
    """将十六进制颜色转换为 RGB"""
    rgb = mcolors.hex2color(hex_color)
    return tuple(int(c * 255) for c in rgb)  # 转换为 0-255 范围
# 示例
hex_color = "#ff5733"  # 示例颜色
rgb_color = hex_to_rgb(hex_color)
print(f"十六进制颜色 {hex_color} 转换为 RGB: {rgb_color}")

valuergb = mcolors.hex2color(hex_color)





## 添加channels
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/pro/

## 配置micromamba/miniconda/channels
![https://mp.weixin.qq.com/s/qt9XvLvHniwQCe6rZAzEZA](https://mp.weixin.qq.com/s/qt9XvLvHniwQCe6rZAzEZA)


## cat ~/.condarc 
auto_activate_base: false
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch-lts: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  deepmodeling: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/



## 颜色 color
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

def plot_colortable(colors, *, ncols=4, sort_colors=True):
    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12
    # Sort colors by hue, saturation, value and name.
    if sort_colors is True:
        names = sorted(
            colors, key=lambda c: tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(c))))
    else:
        names = list(colors)
    n = len(names)
    nrows = math.ceil(n / ncols)
    width = cell_width * ncols + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72
    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height, (width-margin)/width, (height-margin)/height)
    ax.set_xlim(0, cell_width * ncols)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height
        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7
        ax.text(text_pos_x, y, name, fontsize=14, horizontalalignment='left', verticalalignment='center')
        ax.add_patch(Rectangle(xy=(swatch_start_x, y-9), width=swatch_width, height=18, facecolor=colors[name], edgecolor='0.7'))
    return fig
plot_colortable(mcolors.BASE_COLORS, ncols=3, sort_colors=False)
plot_colortable(mcolors.TABLEAU_COLORS, ncols=2, sort_colors=False)
plot_colortable(mcolors.CSS4_COLORS)
plot_colortable(mcolors.XKCD_COLORS)
