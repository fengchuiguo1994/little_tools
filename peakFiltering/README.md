先判断直上直下的区域
```
python fltpeak.py k562.all.ATAC.bw k562.all.ATAC_peaks.narrowPeak
```
然后用igv批量截取这些区域的图像
```
bedtools igv -i k562.all.ATAC_peaks.narrowPeak.drop.txt > test.sh
bedtools igv -i k562.all.ATAC_peaks.narrowPeak.start.txt > test.sh
bedtools igv -i k562.all.ATAC_peaks.narrowPeak.end.txt > test.sh
```
编辑需要保留的区域到文件，再次检查一遍，确定需要被过滤掉的区域。
```
sed 's/_/\t/g' k562.all.ATAC_peaks.narrowPeak.drop.con.txt | awk '{if(NR==FNR){aa[$1"\t"$2"\t"$3]}else{if(($1"\t"$2"\t"$3 in aa)==0){print $0}}}' - k562.all.ATAC_peaks.narrowPeak.drop.txt | bedtools igv -i > test.sh


```