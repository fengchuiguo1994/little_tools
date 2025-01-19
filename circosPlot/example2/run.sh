# grep -v -e hsM -e hsY  ../data/karyotype/karyotype.human.hg38.txt | grep -v band | awk '{$4="chr"$4; print $0}' | sed 's/hs/chr/g' | sed 's/chrx/chrX/g' > karyotype.human.hg38.txt
# bedtools makewindows -g ~/Genome/hg38/hg38.fa.sizes -w 500000 | grep -v -e chrM -e chrY | bedtools intersect -a - -b ~/Genome/hg38/hg38.info -wao | awk '$NF>0' | cut -f 1-3 | uniq -c | awk -v OFS="\t" '{print $2,$3,$4,$1}' > genes_num.txt

# awk '$4>0' /data/home/ruanlab/huangxingyu/CSCC/integrationLoci/HPV.all.CSCC.ex1kb.overlap.bedgraph > HPV.all.CSCC.ex1kb.overlap.bedgraph
# awk '$4>0' /data/home/ruanlab/huangxingyu/CSCC/integrationLoci/HPV.all.CSCC.ex5kb.overlap.bedgraph > HPV.all.CSCC.ex5kb.overlap.bedgraph
