# bedtools bamtobed -bedpe -i ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.singlelinker.paired.UU.nr.bam > test/SCG0192.bedpe
# cp ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.e500.clusters.cis.chiasig.gz .

# python fltbedpe.py SCG0192.bedpe SCG0192.flt.bedpe > SCG0192.flt.bedpe.log

# head -n 7000000 SCG0192.flt.bedpe > SCG0192.7000000.bedpe
for i in 7000 14000 35000 70000 140000 350000 700000 1400000 3500000
do
    # sort -R SCG0192.7000000.bedpe | head -n ${i} > SCG0192.${i}.bedpe
    sleep 0.001
done





samtools view -h ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.singlelinker.paired.UU.nr.bam | head -10024 | samtools view -Sb -o test.bam -
bedtools bamtobed -bedpe -i test.bam > test.bedpe
/data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/cpu-dir/cpu cluster -s 8000 -M -B 1000 -5 5,-20 -3 5,480 -t 32 -O test.e500 -j -x -v 4 -g test.bam 

python bedpe2ipet.py test.bedpe test.bedpe2ipet
sort -k1,1 -k4,4 -k2,2n -k3,3n -k5,5n -k6,6n test.bedpe2ipet > test.sorted.bedpe2ipet
python fltbedpe.py test.sorted.bedpe2ipet test.sorted.flt.bedpe2ipet > test.sorted.flt.bedpe2ipet.log






/data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/cpu-dir/cpu cluster -s 8000 -M -B 1000 -5 5,-20 -3 5,480 -t 2 -O test1.e500 -j -x -v 4 -g ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.singlelinker.paired.UU.nr.bam > test1.log 2>&1
/data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/cpu-dir/cpu cluster -s 8000 -m -B 1000 -5 5,0 -3 3,500   -t 2 -O test2.e500 -j -x -v 4 -g ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.singlelinker.paired.UU.nr.bam > test2.log 2>&1

/data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/cpu-dir/cpu cluster -s 8000 -M -B 1000 -5 5,-250 -3 5,250 -t 2 -O test3.e500 -j -x -v 4 -g ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.singlelinker.paired.UU.nr.bam > test3.log 2>&1


/data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/cpu-dir/cpu cluster -s 8000 -M -B 1000 -5 5,-250 -3 5,50 -t 2 -O test4.e500 -j -x -v 4 -g ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.singlelinker.paired.UU.nr.bam > test4.log 2>&1

/data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/cpu-dir/cpu-dir/cpu cluster -s 8000 -M -B 1000 -5 5,-250 -3 5,50 -t 2 -O test5.e500 -j -x -v 4 -g ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.singlelinker.paired.UU.nr.bam > test5.log 2>&1
/data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/cpu-dir/cpu span -s 8000 -t 2 ../SCG0192_GT22-15872_SI-NA-D6_bulk/SCG0192_GT22-15872_SI-NA-D6_bulk.singlelinker.paired.UU.nr.bam > span.log 2>&1
