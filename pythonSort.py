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