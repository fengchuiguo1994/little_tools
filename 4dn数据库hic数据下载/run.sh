# singularity run /data/home/ruanlab/huangxingyu/Tools/hicexplorer.sif hicConvertFormat -m 4DNFITUOMFUQ.hic -o 4DNFITUOMFUQ.mcool --inputFormat hic --outputFormat cool
# zcat GSE63525_K562_HiCCUPS_looplist.txt.gz | sed '1d' | awk -v OFS="\t" '{print "chr"$1,$2,$3,"chr"$4,$5,$6,$8}' > GSE63525_K562_HiCCUPS.bedpe

# cat <(zcat 4DNFI2R1W3YW.pairs.gz | head -100000 | grep "^#") <(zcat 4DNFI2R1W3YW.pairs.gz | grep -v "^#" | awk '$8=="UU" && $9+0.0>=20 && $10+0.0>=20') > 4DNFI2R1W3YW.flt.pairs # 907137059 # 726879957
# sbatch -p ruan_cpu -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=100G --wrap=" singularity run /data/home/ruanlab/huangxingyu/scChAIR/ChAIR.sif juicer_tools pre -r 2500000,2000000,1000000,500000,250000,200000,100000,50000,25000,20000,10000,5000,2500,2000,1000 4DNFI2R1W3YW.flt.pairs 4DNFI2R1W3YW.flt.pairs.hic hg38 "

# cat 4DNFI2R1W3YW.flt.pairs | grep -v "^#" | awk -v OFS="\t" '{if($4=="+"){start1=$3; end1=$3+100;}else{start1=$3-100; end1=$3;} if($7=="+"){start2=$5; end2=$5+100;}else{start2=$5-100; end2=$5;} mapq=$9; if($9>$10){mapq=$10;} print $2,start1,end1,$4,start2,end2,$1,mapq,$6,$7}' > 4DNFI2R1W3YW.bedpe




# awk -v OFS="\t" '{print $1,$2,$3,$7,$8,$9"\n"$4,$5,$6,$7,$8,$10}' 4DNFI2R1W3YW.bedpe | bedtools bedtobam -i - -g ~/Genome/hg38/hg38.fa.sizes > 4DNFI2R1W3YW.bam
# sbatch -p ruan_cpu -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=20G --wrap=" samtools sort -@ 10 -o 4DNFI2R1W3YW.sorted.bam 4DNFI2R1W3YW.bam "


# sbatch -p ruan_cpu -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 2 --mem=200G --wrap=" python bedpe2bed.py 4DNFI2R1W3YW.bedpe 4DNFI2R1W3YW.bedpe.tobed "
# for i in {1..22} X Y; do sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=20G --wrap=" bedtools bedtobam -i 4DNFI2R1W3YW.bedpe.tobed.chr${i}.bed -g ~/Genome/hg38/hg38.fa.sizes | samtools sort -@ 10 -o 4DNFI2R1W3YW.bedpe.tobed.chr${i}.bam - "; done
# samtools merge -@ 10 -o 4DNFI2R1W3YW.bedpe.tobed.bam $(ls 4DNFI2R1W3YW.bedpe.tobed.chr*.bam)

# sbatch -p ruan_cpu -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=20G --wrap=" samtools index 4DNFI2R1W3YW.bedpe.tobed.bam && bamCoverage -b 4DNFI2R1W3YW.bedpe.tobed.bam -o 4DNFI2R1W3YW.bedpe.tobed.bw -p 10 -bs 5 "

