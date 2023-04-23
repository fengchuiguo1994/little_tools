bash star.align.pair.sh > star.align.pair.sh.list
awk '{print $1"\t"$NF}' star.align.pair.sh.list | while read a b; do seff $b > $a.run.sh.log; done