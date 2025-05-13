## bwa+pairtools
```
pairtools merge -o cancer.pairs.gz CH001.dedup.pairs.gz.pairs.cut.gz CH005.dedup.pairs.gz.pairs.cut.gz CH008.dedup.pairs.gz.pairs.cut.gz CH009.dedup.pairs.gz.pairs.cut.gz CH011.dedup.pairs.gz.pairs.cut.gz CH013.dedup.pairs.gz.pairs.cut.gz CH015.dedup.pairs.gz.pairs.cut.gz CH017.dedup.pairs.gz.pairs.cut.gz CH019.dedup.pairs.gz.pairs.cut.gz
java -Xmx48000m -Djava.awt.headless=true -jar /data/home/ruanlab/huangxingyu/Tools/juicer_tools.1.8.9_jcuda.0.8.jar pre -r 2500000,2000000,1000000,500000,250000,200000,100000,50000,25000,20000,10000,5000,2500,2000,1000 cancer.pairs.gz cancer.pairs.hic hg38
/data/home/ruanlab/huangjiaxiang/miniconda3/envs/hicexplorer/bin/hicConvertFormat -m cancer.nochr.pairs.hic --inputFormat hic --outputFormat cool -o cancer.nochr.pairs.hic2cool
for i in `cooler ls cancer.nochr.pairs.hic2cool.mcool`;  do cooler balance --max-iters 500 $i; done
cooltools expected-cis --smooth --aggregate-smoothed cancer.nochr.pairs.hic2cool.mcool::/resolutions/1000 --nproc 5 -o cancer.nochr.pairs.hic2cool.mcool.1k.tsv
```

## corrlation
```
python hicrep.auto.py  -i CH001.dedup.pairs.gz.pairs.cut.hic.mcool,CH002.dedup.pairs.gz.pairs.cut.hic.mcool,CH005.dedup.pairs.gz.pairs.cut.hic.mcool -o CSCC.hic.cor -s hg38.hic.size
Rscript pheatmap.r CSCC.hic.cor.cor.mat CSCC.hic.cor.cor.mat.pdf CSCC.hic.cor.cor.mat.cluster.pdf
```

## loop+TAD+compartment
```
for i in {1..22} X; do java -jar ~/Tools/juicer_tools.1.8.9_jcuda.0.8.jar eigenvector -p KR normal.pairs.hic chr${i} BP 1000000 normal.chr${i}.abcompartment; java -jar ~/Tools/juicer_tools.1.8.9_jcuda.0.8.jar eigenvector -p KR cancer.pairs.hic chr${i} BP 1000000 cancer.chr${i}.abcompartment; done
grep -v -e "chrY" -e "chrM" /data/home/ruanlab/huangxingyu/Genome/hg38/hg38.fa.sizes | bedtools makewindows -g - -w 1000000 -s 1000000 > hg38.1M.bed
# for i in {1..22} X; do cat normal.chr${i}.abcompartment; done | paste hg38.1M.bed - > normal.abcompartment


cooler dump --header -t bins CH001.dedup.pairs.gz.pairs.cut.hic.mcool::resolutions/100000 | cut -f1-3 > bins.100000.tsv
# cooltools genome gc bins.100000.tsv /data/home/ruanlab/huangxingyu/Genome/hg38/hg38.fa > gc.100000.tsv
# cat <(head -1 bins.100000.tsv) <(sed '1d' bins.100000.tsv | awk '{print "chr"$0}') > bins.100000.chr.tsv
# cooltools genome gc bins.100000.chr.tsv /data/home/ruanlab/huangxingyu/Genome/hg38//hg38.fa > gc.100000.tsv
# cooltools eigs-cis -o cancer.eigs.100000 --view cancer.tsv --phasing-track gc.100000.tsv --n-eigs 1 cancer.hic2cool.mcool::resolutions/100000


mkdir normalTAD cancerTAD
# java -jar ~/Tools/juicer_tools.1.8.9_jcuda.0.8.jar arrowhead -r 25000 -k KR normal.pairs.hic normalTAD


hicFindTADs -m cancer.hic2cool.KR.cool --outPrefix cancerTAD/cancerTAD.hicFindTADs.result -p 1 --correctForMultipleTesting fdr


java -jar ~/hicFileForCompareHiCChIA/K562/juicer_tools_1.13.02.jar hiccups --cpu -m 1024 -r 5000,10000 -k KR normal.pairs.hic normal.KR.cpu.loop
```