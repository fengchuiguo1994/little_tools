#!/bin/bash

threads=10
index=/data/home/ruanlab/huangxingyu/Tools/cellranger-arc-2.0.2/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star
species=human
outdir="./"
indir="../"

picard=/data/home/ruanlab/huangxingyu/Tools/littletools/picard.jar

for i in ENCSR000AEM_Rep1 ENCSR000AEM_Rep2 ENCSR000AEO_Rep1 ENCSR000AEO_Rep2 ENCSR384ZXD_Rep1 ENCSR384ZXD_Rep2 ENCSR530NHO_Rep1 ENCSR530NHO_Rep2 ENCSR545DKY_Rep1 ENCSR545DKY_Rep2 ENCSR594NJP_Rep1 ENCSR594NJP_Rep2 ENCSR596ACL_Rep1 ENCSR596ACL_Rep2
do
cat << ENDHERE > ${i}.align.sh
#!/usr/bin/bash

module load R-3.6.1

/data/home/ruanlab/huangxingyu/Tools/cellranger-arc-2.0.1/lib/bin/STAR --runMode alignReads --genomeDir ${index} --readFilesIn ${indir}/${i}_R1.clean.pair.fq.gz ${indir}/${i}_R2.clean.pair.fq.gz --readFilesCommand zcat --runThreadN ${threads} --outFileNamePrefix ${outdir}/${i} --outFilterMultimapNmax 10 --outSAMattributes All --outSAMtype BAM Unsorted --outSAMunmapped Within --limitOutSJcollapsed 32000000 --limitIObufferSize 1200000000 2000000000 --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5
samtools sort -@ ${threads} -o ${outdir}/$i.bam ${outdir}/${i}Aligned.out.bam && samtools index ${outdir}/$i.bam
samtools flagstat ${outdir}/$i.bam > ${outdir}/$i.bam.stat

java -jar ${picard} CollectInsertSizeMetrics -I ${outdir}/${i}.bam -O ${outdir}/${i}.insert_size_metrics.txt -H ${outdir}/${i}.insert_size_histogram.pdf -M 0.5 --VALIDATION_STRINGENCY SILENT 

samtools view -q 30 -o ${outdir}/${i}.flt.bam ${outdir}/${i}.bam
samtools flagstat ${outdir}/${i}.flt.bam > ${outdir}/${i}.flt.bam.stat
samtools index ${outdir}/${i}.flt.bam

bamCoverage -o ${outdir}/$i.RPKM.bw --binSize 5 -b ${outdir}/$i.flt.bam -p ${threads} --normalizeUsing RPKM 
bamCoverage -o ${outdir}/$i.bw --binSize 5 -b ${outdir}/$i.flt.bam -p ${threads}
bamCoverage -o ${outdir}/$i.RPKM.bedgraph --binSize 5 -b ${outdir}/$i.flt.bam -p ${threads} --normalizeUsing RPKM -of bedgraph
bamCoverage -o ${outdir}/$i.bedgraph --binSize 5 -b ${outdir}/$i.flt.bam -p ${threads} -of bedgraph

bamCoverage -o ${outdir}/$i.fwd.bw --binSize 5 -b ${outdir}/$i.flt.bam -p ${threads} --filterRNAstrand forward
bamCoverage -o ${outdir}/$i.rev.bw --binSize 5 -b ${outdir}/$i.flt.bam -p ${threads} --filterRNAstrand reverse
bamCoverage -o ${outdir}/$i.fwd.bedgraph --binSize 5 -b ${outdir}/$i.flt.bam -p ${threads} --filterRNAstrand forward -of bedgraph
bamCoverage -o ${outdir}/$i.rev.bedgraph --binSize 5 -b ${outdir}/$i.flt.bam -p ${threads} --filterRNAstrand reverse -of bedgraph

cat <(awk -v OFS="\t" '{if($4!=0){print $1,$2,$3,$4,0}}' ${outdir}/$i.fwd.bedgraph) <(awk -v OFS="\t" '{if($4!=0){print $1,$2,$3,0,$4}}' ${outdir}/$i.rev.bedgraph) > ${outdir}/$i.forBASIC.ssbedgraph
ENDHERE

rid=$(sbatch -J ${i} --mem=40G -e ${i}.err -o ${i}.out -N 1 -n ${threads} -p cpu ${i}.align.sh)
# echo -e "${i}\t${rid}"
echo -e "${i}\t${rid}" | awk '{print $1"\t"$NF}'
done

