## 多线程测试
import multiprocessing
from multiprocessing import Manager
import time
import random
import copy

mydata = Manager().list()

def func(msg,nn):
    print(msg)
    print(nn)
    newlist = copy.deepcopy(mydata)
    random.shuffle(newlist)
    time.sleep(1)
    print(newlist)
    return "123"
    
if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=4)
    # mydata = [1,2,3,4,5,6,7,8,9,10]
    mydata.append(1)
    mydata.append(2)
    mydata.append(3)
    los = ((1,2),(2,3),(3,4),(4,5),(5,6),(6,7))
    results = pool.map(func,los)
    
    print("mark")
    pool.close()
    pool.join()
    print("finish")
    for result in results:
        print(result)



## 单行命令
#### 对于gencode(类似)的注释文件快速提取gene区间和基因的信息
```
perl -lane 'if($F[2]eq"gene"){if(/gene_id "(.+?)\.\d+";.+?gene_type "(.+?)";.+?gene_name "(.+?)";/){print "$F[0]\t".($F[3]-1)."\t$F[4]\t$1\t.\t$F[6]\t$3\t$2"}}' gencode.vM23.annotation.gtf > mm10.info
```