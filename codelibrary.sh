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