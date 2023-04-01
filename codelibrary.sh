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