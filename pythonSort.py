# 列表排序
out = []
out.sort(key=lambda x:(x.chrom,x.start,x.end)) # 对象元素排序
out.sort(key=lambda x: (x[0],x[1])) # 列表元素排序
out.sort(key=lambda x:(x[0].chrom,x[1].chrom,x[0].start,x[1].start,x[0].end,x[1].end)) # 列表元素，列表元素中含有对象元素，排序


import re
out.sort(key=lambda x:(int(re.search("\d+",x[0])[0]),x[1])) # 列表元素，列表元素中的字符串元素中的数字部分，排序


# 数字和字符串混和排序，在python3中取消了内建函数cmp，提供了functools来转化
import functools
def cmp(a, b):
    if b < a:
        return 1
    if a < b:
        return -1
    return 0
def mycmp(o1,o2):
    if isinstance(o1,str) and isinstance(o2,str): 
        return cmp(o1,o2)
    elif isinstance(o1,int) and isinstance(o2,int): 
        return cmp(o1,o2)
    elif isinstance(o1,str) and isinstance(o2,int): 
        return 1
    else:
        return -1
a = [1, 2, 5, 4]
print(sorted(a, key=functools.cmp_to_key(cmp)))
aa = [5,4,"5","4",3,"11",2,"3","2",1,"1"]
print(aa)
aa.sort(key=functools.cmp_to_key(mycmp))
print(aa)


import functools
def cmp(a, b):
    if b < a:
        return 1
    if a < b:
        return -1
    return 0
def mycmp(o1,o2):
    o1 = o1.replace("chr","")
    o2 = o2.replace("chr","")
    if o1.isdigit() and o2.isdigit():
        return cmp(int(o1),int(o2))
    elif not o1.isdigit() and o2.isdigit():
        return 1
    elif o1.isdigit() and not o2.isdigit():
        return -1
    else:
        return cmp(o1,o2)
aa = ["chrM","chr1","chrY","HPV18","chrX","chr3","chr2"]
# aa = ["chr1","chr12","chr3","chr2","chr11"]
aa.sort(key=functools.cmp_to_key(mycmp))
print(aa)


import functools
import re
def cmp(a, b):
    if b < a:
        return 1
    if a < b:
        return -1
    return 0
def blockcmp(o1,o2):
    reo1 = re.search("([-a-zA-Z]+)(\d+)",o1)
    reo2 = re.search("([-a-zA-Z]+)(\d+)",o2)
    if reo1[1] == reo2[1]:
        return cmp(int(reo1[2]),int(reo2[2]))
    else:
        return cmp(reo1[1],reo2[1])
    
aa = "WGS-block30057,WGS-block7729,WGS-block28602,HIC-block5052,WGS-block39719,WGS-block34376,WGS-block13793,WGS-block9085,WGS-block2935,WGS-block37164,WGS-block39540,WGS-block16211,WGS-block11248,WGS-block37433,WGS-block40022,WGS-block5263,WGS-block1790,WGS-block34214,WGS-block34671,WGS-block4724,WGS-block9860,WGS-block7994,WGS-block23228,WGS-block3813,WGS-block13748,WGS-block12850,WGS-block21163,WGS-block25924,WGS-block26696,WGS-block12219,WGS-block25865,WGS-block35742,WGS-block40630,WGS-block9667,WGS-block36017,WGS-block22469,WGS-block31982,WGS-block17917,WGS-block351,HIC-block2023,WGS-block40653,WGS-block14945,HIC-block640,WGS-block17087,WGS-block1984,HIC-block2489,WGS-block38988,WGS-block12976,WGS-block12404,WGS-block18985,WGS-block38831,WGS-block25118,WGS-block40535,WGS-block39288,WGS-block40333,WGS-block19381,HIC-block5505,WGS-block27006,WGS-block40261,WGS-block23895,WGS-block38273,WGS-block1593,WGS-block25128,WGS-block8895,WGS-block26330,WGS-block34551,WGS-block23914,WGS-block5692,WGS-block3592,WGS-block19126,WGS-block26263,WGS-block29872,WGS-block15680,WGS-block29716,WGS-block16692,WGS-block19995,WGS-block27259,WGS-block30760,WGS-block31619,WGS-block36218,WGS-block22722,WGS-block17310,WGS-block5381,WGS-block9577,WGS-block12010,WGS-block35898,WGS-block10752,WGS-block12260,WGS-block38777,WGS-block18490,WGS-block39621,WGS-block8090,WGS-block3481,WGS-block36598,WGS-block8410,WGS-block40227,HIC-block1797"
bb = aa.split(",")
bb.sort(key=functools.cmp_to_key(blockcmp))
for i in bb:
    print(i)


## 字典排序（转化成列表再排序的）
aa = {1:5,2:4,3:3,4:2,5:1} 
aabb = sorted(aa.items(),key=lambda x:x[0]) # 按照键值对排序
aabb = sorted(aa.items(),key=lambda x:x[1]) # 按照键值对排序