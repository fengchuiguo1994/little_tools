## shell控制的循环
案例一：
```
total=100
nthread=7

let per=$total/$nthread
echo $per
let rest=$total-$per*$nthread
echo $rest

index=1
n=0
while ((index<=nthread))
do
    if ((index<=rest))
    then
        echo -n "$n"
        let n=$n+per+1
        echo " $n"
    else
        echo -n "$n"
        let n=$n+per
        echo " $n"
    fi
    let index=$index+1
done
```
案例二
```
# testrun.sh  
echo $1
sleep 10
echo $1

# testtest.sh
NTHREADS=2
function jobmax
{
    typeset -i MAXJOBS=${NTHREADS}
    sleep .2
    while (( ($(pgrep -P $$ | wc -l) - 1) >= $MAXJOBS ))
    do
        sleep 5
    done
}

for i in {1..20}
do
    bash testrun.sh $i &
    jobmax
done
```

## 帮助文档
SOFT="loaddata.sh"
VERSION="1.0.0"

function usage {
    echo -e "usage: bash loadbedpe.sh -i inputfile -g genome -f uploadfolder -t datatype"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "$SOFT $VERSION"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT : inputfile"
    echo "   -g|--genome GENOME : Must be a genome that exists on the BASIC browser"
    echo "   -f|--folder FOLDER : Upload to the specified library"
    echo "   -t|--datatype DATATYPE : Type of data to be uploaded [bedgraph, bed, ]"
    echo "   [-v|--version]: version"
    echo "   [-h|--help]: help"
    exit;
}


function version {
    echo -e "$SOFT version $VERSION"
    exit
}

if [ $# -lt 1 ]
then
    usage
    exit
fi

while [ -n "$1" ]
do 
    case $1 in
        -i|--inputfile) 
            input=$2
            shift 2
            ;;
        -g|--genome)
            genome=$2
            shift 2
            ;;
        -f|--folder)
            folder=$2
            shift 2
            ;;
        -t|--datatype)
            datatype=$2
            shift 2
            ;;
        -h|--help)
            help
            exit 0
            ;;
        -v|--version)
            version
            exit 0
            ;;
        --) 
            shift
            break
            ;;
        -*) 
            echo "error: no such option $1."
            exit 1 
            ;;
        # *) echo "error: no such option $1."; exit 1;;
        *) 
            break
            ;;
     esac
done

echo $1
echo $2

echo $input
echo $genome
echo $folder
echo $datatype

## 递归遍历得到所有文件
#!/bin/bash

# 采集一个函数
readDir() {
  # 获取传入的目录路径
  local dir=$1
  # 循环指定目录下的所有文件
  local files
  files=$(ls "$dir")
  for file in $files; do
    local path="$dir/$file" #指的是当前遍历文件的完整路径
    # 判断是否是目录，如果是目录则递归遍历，如果是文件则打印该文件的完整路径
    if [ -d "$path" ]; then
      readDir "$path"
    else
      echo "$path"
    fi
  done
}

readDir root

## 调用系统参数
uname -a

cat /etc/*-release | sort -u

cat /proc/version

glibc version
ldd --version | head -n 1

grep -m 1 'model name' /proc/cpuinfo | cut -d ':' -f 2 | sed 's/^[ \t]*//'

grep 'physical id' /proc/cpuinfo | sort -u | wc -l

grep -c processor /proc/cpuinfo

grep -m 1 'flags' /proc/cpuinfo | cut -d ':' -f 2 | sed 's/^\s*//'

grep MemTotal /proc/meminfo | cut -d ':' -f 2 | sed 's/^[ \t]*//'

df -Ph | awk '{print $2, $3, $4}'

mount | cut -d ' ' -f 5,6

bash -c 'ulimit -a'

bash -c 'ulimit -aH'

cat /proc/sys/fs/file-{max,nr}

sysctl vm

cat /sys/kernel/mm/*transparent_hugepage/enabled

cat /proc/self/cgroup

cgroup mem stats

cat /sys/fs/cgroup/memory/user.slice/memory.*soft_limit_in_bytes

cat /sys/fs/cgroup/memory/user.slice/memory.limit_in_bytes

cat /sys/fs/cgroup/memory/user.slice/memory.memsw.limit_in_bytes

head -n 1 /proc/1/sched | cut -d ' ' -f 1

which qsub

which bsub

which configureBclToFastq.pl

which bcl2fastq

which java

java -version 2>&1 | cat

echo $TENX_REFDATA

cat $TENX_REFDATA/version

which qconf

sinfo -O nodes,maxcpuspernode,memory,time

mrp --version

ls $(dirname $(dirname $(which mrp)))/jobmanagers/*.template

## 计算程序运行时间和内存消耗
```
time
/usr/bin/time

/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' grep -f file.txt pattern.txt
# 注意这里的内存单位为kb
```

## ImageMagick图片处理
格式转换
```
convert foo.jpg foo.png
find ./ -name "*.jpg"|xargs -i basename {} .jpg|xargs -i convert {}.jpg {}.png

convert -resize 50%x50%  -quality 70 -strip foo.jpg bar.jpg
find ./ -name "*.jpg"|xargs -i basename {} .jpg|xargs -i convert -resize 50%x50% {}.jpg {}_50.jpg

convert -blur 80 foo.jpg foo.png # 高斯模糊
convert -monochrome foo.png bar.png # 把图片变为黑白颜色
convert -paint 4 foo.png bar.png # 油画效果
convert -charcoal 2 foo.png bar.png # 铅笔画效果
```

## 安装包
#### 安装gcc
[gcc需要的三个gmp、mpc、mpfr包](https://gcc.gnu.org/pub/gcc/infrastructure/)，[下载链接2](https://ftp.gnu.org/gnu/)
```
cd Tools
make gcc
make gcc/gmp gcc/mpc gcc/mpfr

tar -jxvf gmp-4.3.2.tar.bz2
cd gmp-4.3.2
./configure --prefix=/home/Tools/gcc/gmp/ #gmp安装路径
make
make install
cd ..

tar -jxvf mpfr-2.4.2.tar.bz2
cd mpfr-2.4.2
./configure --prefix=/home/Tools/gcc/mpfr/ --with-gmp=/home/Tools/gcc/gmp/ #congfigure后面是mpfr安装路径及依赖的gmp路径
make
make install
cd ..

tar -zxvf mpc-0.8.1.tar.gz
cd mpc-0.8.1
./configure --prefix=/home/Tools/gcc/mpc/ --with-gmp=/home/Tools/gcc/gmp/ --with-mpfr=/home/Tools/gcc/mpfr/
make
make install
cd ..

添加环境变量到~/.bashrc
export LD_LIBRARY_PATH=/home/Tools/gcc/gmp/lib/:/home/Tools/gcc/mpfr/lib/:/home/Tools/gcc/mpc/lib/
export LIBRARY_PATH=$LD_LIBRARY_PATH
```
[安装gcc](https://gcc.gnu.org/releases.html)
```
make gcc/gcc
tar -jxvf gcc-8.4.0.tar.bz2
cd gcc-gcc-8.4.0
./configure --prefix=/home/Tools/gcc/gcc/ --enable-threads=posix --disable-checking --disable-multilib --with-mpc=/home/Tools/gcc/mpc/ --with-gmp=/home/Tools/gcc/gmp/ --with-mpfr=/home/Tools/gcc/mpfr/ 
make -j 10
make install
```