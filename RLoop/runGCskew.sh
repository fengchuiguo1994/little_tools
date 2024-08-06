#perl -lane 'if(/>/){$_=~s/ .+?$//;print $_}else{print $_}' MH63RS2.LNNK00000000.fasta > MH63RS2.fa
bedtools makewindows -g MH63RS2.size -w 200 -s 10 > bins.bed
bedtools nuc -fi MH63RS2.LNNK00000000.fasta -bed bins.bed > bins.count.bed
grep -v "#" bins.count.bed | awk '{cg=0;at=0;}{if($8+$7>0){cg=($8-$7)/($8+$7);}if($6+$9>0){at=($6-$9)/($6+$9);}print $1"\t"$2"\t"$3"\t"cg"\t"at;}' > GCATskew.bed 
awk '{if($3-$2==200){ee=int(($3+$2)/2);print $1"\t"ee-5"\t"ee+5"\t"$4}}' GCATskew.bed > GC.bedgraph
~/software/bedGraphToBigWig GC.bedgraph MH63RS2.size GC.bw
awk '{if($3-$2==200){ee=int(($3+$2)/2);print $1"\t"ee-5"\t"ee+5"\t"$5}}' GCATskew.bed > AT.bedgraph
~/software/bedGraphToBigWig AT.bedgraph MH63RS2.size AT.bw
computeMatrix scale-regions -S AT.bw GC.bw -R MH63RS2.bed -a 2000 -b 2000 -m 2000 --skipZeros --outFileName gene.richheatmap.gz --numberOfProcessors 4
zcat gene.richheatmap.gz | perl -lane 'if($.==1){print $_}else{if($F[5] eq "+"){print $_}else{foreach(6..$#F){$F[$_]=-$F[$_];}print join("\t",@F)}}' | gzip > gene.richheatmap.tmp.gz
plotProfile -m gene.richheatmap.tmp.gz -out gene.richheatmap.tmp.png --perGroup
