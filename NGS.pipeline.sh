#### RNA-Seq
## Fastqc、Trimmomatic、hisat2、RSeQC、samtools、htseq、picard、deeptools
thread=10
prefix=out
outdir=./
knownspl=genome/genome.ss # python hisat2dir/extract_splice_sites.py genome.gtf > genome.ss
index=genome/genome.hisat # python hisat2dir/extract_exons.py genome.gtf > genome.exon && hisat2-build -p 8 --ss genome.ss --exon genome.exon genome.fa genome.hisat
species=human
gtffile=gencode.V19.gtf # null
gfffile=gencode.V19.gff # null
exon_gene_length=hg38.exon_gene_length # perl -F"\t" -lane 'if($F[2] eq "exon"){$gene="tmp";$gene=$1 if($F[8]=~/gene_id "(.+?)";/);print "$gene#$F[6]#$F[0]\t".($F[3]-1)."\t$F[4]"}' gencode.V19.gtf | bedtools sort -i - | bedtools merge -i - | perl -F"\t" -lane '@tmp=split(/#/,$F[0]);print "$tmp[2]\t$F[1]\t$F[2]\t$tmp[1]\t.\t$tmp[0]"' > gencode.V19.exon_gene_length
bedfile=gencode.V19.bed12
IN1=fq1
IN2=fq2


Trimmomaticdir=/public/home/xyhuang/Tools/Trimmomatic-0.36/
Trimmomaticjar=trimmomatic-0.36.jar
Picarddir=/public/home/xyhuang/Tools/javatools/
Tooldir=/public/home/xyhuang/Tools/littletools

mkdir -p ${outdir}/rawqc
fastqc -o ${outdir}/rawqc -t ${thread} ${IN1} ${IN2}

mkdir -p ${outdir}/trim
java -jar ${Trimmomaticdir}/${Trimmomaticjar} PE -threads ${thread} -phred33 ${IN1} ${IN2} ${outdir}/trim/${prefix}_R1.fq.gz ${outdir}/trim/${prefix}_R1.unpair.fq.gz ${outdir}/trim/${prefix}_R2.fq.gz ${outdir}/trim/${prefix}_R2.unpair.fq.gz ILLUMINACLIP:${Trimmomaticdir}/adapters/TruSeq3-PE-2.fa:2:30:7:8:true LEADING:10 TRAILING:10 AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:30 && fastqc -t ${thread} ${outdir}/trim/${prefix}_R1.fq.gz ${outdir}/trim/${prefix}_R1.unpair.fq.gz ${outdir}/trim/${prefix}_R2.fq.gz ${outdir}/trim/${prefix}_R2.unpair.fq.gz

# 需要首先判断是否是连特异性以及方向，但一般适用于RF
# hisat2 -x ${index} -1 ${outdir}/trim/${prefix}_R1.fq.gz -2 ${outdir}/trim/${prefix}_R2.fq.gz | head -2000000 > ${outdir}/align/${prefix}.tmp.test.sam.tmp
infer_experiment.py -i ${outdir}/align/${prefix}.tmp.test.sam.tmp -r ${bedfile} -s 1000000 > ${outdir}align/${prefix}.sam.strand &

mkdir -p ${outdir}/align
hisat2 --known-splicesite-infile ${knownspl} --rna-strandness RF --dta --rg-id ${species} --rg LB:totalRNA --rg PG:hisat2  --rg PL:ILLUMINA --rg PU:lane --rg SM:${prefix} -p ${thread} -x ${index} -1 ${outdir}/trim/${prefix}_R1.fq.gz -2 ${outdir}/trim/${prefix}_R2.fq.gz | tee ${outdir}align/${prefix}.sam | samtools view -Sb - | samtools sort -@ ${thread} -o ${outdir}/align/${prefix}.bam -

# python ${Tooldir}/SplitHisatAlignFile.py ${outdir}align/${prefix}.bam ${outdir}align/${prefix} && cp ${outdir}align/${prefix}.uniq-uniq.bam ${outdir}align/${prefix}.flt.bam
samtools view -Sb -q 20 -o ${outdir}align/${prefix}.flt.bam ${outdir}align/${prefix}.bam
java -jar ${Picarddir}/picard.jar CollectInsertSizeMetrics I=${outdir}align/${prefix}.flt.bam O=${outdir}align/${prefix}.insert_size_metrics.txt H=${outdir}align/${prefix}.insert_size_histogram.pdf M=0.5 VALIDATION_STRINGENCY=SILENT
samtools index ${outdir}align/${prefix}.flt.bam && samtools flagstat ${outdir}align/${prefix}.bam > ${outdir}align/${prefix}.bam.stat && samtools flagstat ${outdir}align/${prefix}.flt.bam > ${outdir}align/${prefix}.flt.bam.stat

bamCoverage -o ${outdir}align/${prefix}.RPKM.bw -bs 5 -b ${outdir}align/${prefix}.flt.bam -p ${thread} --minMappingQuality 20 --normalizeUsing RPKM && bamCoverage -o ${outdir}align/${prefix}.bw -bs 5 -b ${outdir}align/${prefix}.flt.bam -p ${thread} --minMappingQuality 20
bamCoverage -o ${outdir}align/${prefix}.fwd.RPKM.bw --filterRNAstrand forward -bs 5 -b ${outdir}align/${prefix}.flt.bam -p ${thread} --minMappingQuality 20 --normalizeUsing RPKM && bamCoverage -o ${outdir}align/${prefix}.fwd.bw --filterRNAstrand forward -bs 5 -b ${outdir}align/${prefix}.flt.bam -p ${thread} --minMappingQuality 20
bamCoverage -o ${outdir}align/${prefix}.rev.RPKM.bw --filterRNAstrand reverse -bs 5 -b ${outdir}align/${prefix}.flt.bam -p ${thread} --minMappingQuality 20 --normalizeUsing RPKM && bamCoverage -o ${outdir}align/${prefix}.rev.bw --filterRNAstrand reverse -bs 5 -b ${outdir}align/${prefix}.flt.bam -p ${thread} --minMappingQuality 20

if [ ${gtffile} != 'null' ];then
    htseq-count -f sam -r name -s reverse -t gene -i gene_id -m intersection-nonempty ${outdir}align/${prefix}.sam ${gtffile} > ${outdir}align/${prefix}.gtf.gene.count &
    htseq-count -f sam -r name -s reverse -t exon -i gene_id -m intersection-nonempty ${outdir}align/${prefix}.sam ${gtffile} > ${outdir}align/${prefix}.gtf.exon.count &
    wait
    python ${Tooldir}/count2FPKM.py ${outdir}align/${prefix}.gtf.gene.count ${exon_gene_length} ${outdir}align/${prefix}.gtf.gene.count.fpkm
fi
if [ ${gfffile} != 'null' ];then
    htseq-count -f sam -r name -s reverse -t gene -i ID -m intersection-nonempty ${outdir}align/${prefix}.sam ${gtffile} > ${outdir}align/${prefix}.gff.gene.count &
    htseq-count -f sam -r name -s reverse -t exon -i ID -m intersection-nonempty ${outdir}align/${prefix}.sam ${gtffile} > ${outdir}align/${prefix}.gff.exon.count &
    wait
    python ${Tooldir}/count2FPKM.py ${outdir}align/${prefix}.gff.gene.count ${exon_gene_length} ${outdir}align/${prefix}.gff.gene.count.fpkm
fi


# multiqc -o ${outdir}/all_combine_raw_multi ${outdir}/rawqc
# multiqc -o ${outdir}/all_combine_filter.multi ${outdir}/trim



#### ChIP-Seq
## Fastqc、Trimmomatic、bwa、samtools、picard、deeptools
thread=10
prefix=out
outdir=./
index=genome.fa # bwa index genome.fa
GENOME_LEN=36000000 # perl -lane '{$_=~s/[Nn]+//g;$sum+=length($_);}END{print $sum;}' hg38.fa 未考虑repeat的影响，但是如果所有的都设置一样的话就没有影响了。
species=human
IN1=fq1
IN2=fq2
input.bam=inputbam

Trimmomaticdir=/public/home/xyhuang/Tools/Trimmomatic-0.36/
Trimmomaticjar=trimmomatic-0.36.jar
Picarddir=/public/home/xyhuang/Tools/javatools/
Tooldir=/public/home/xyhuang/Tools/littletools
ppddir=/public/home/xyhuang/Tools/phantompeakqualtools/

mkdir -p ${outdir}/rawqc
fastqc -o ${outdir}/rawqc -t ${thread} ${IN1} ${IN2}

mkdir -p ${outdir}/trim
java -jar ${Trimmomaticdir}/${Trimmomaticjar} PE -threads ${thread} -phred33 ${IN1} ${IN2} ${outdir}/trim/${prefix}_R1.fq.gz ${outdir}/trim/${prefix}_R1.unpair.fq.gz ${outdir}/trim/${prefix}_R2.fq.gz ${outdir}/trim/${prefix}_R2.unpair.fq.gz ILLUMINACLIP:${Trimmomaticdir}/adapters/TruSeq3-PE-2.fa:2:30:7:8:true LEADING:10 TRAILING:10 AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:30 && fastqc -t ${thread} ${outdir}/trim/${prefix}_R1.fq.gz ${outdir}/trim/${prefix}_R1.unpair.fq.gz ${outdir}/trim/${prefix}_R2.fq.gz ${outdir}/trim/${prefix}_R2.unpair.fq.gz

mkdir -p ${outdir}/align
bwa mem -t ${thread} $index ${outdir}/trim/${i}_R1.fq.gz ${outdir}/trim/${i}_R2.fq.gz -R "@RG\tID:${prefix}\tSM:${prefix}\tLB:${species}\tPL:illumina\tPU:run" | samtools view -Sb - | samtools sort -@ ${thread} -o ${outdir}/align/${prefix}.bam - && samtools index ${outdir}/align/${prefix}.bam && samtools flagstat ${outdir}/align/${prefix}.bam > ${outdir}/align/${prefix}.bam.stat

java -jar ${Picarddir}/picard.jar CollectInsertSizeMetrics I=${outdir}/align/${prefix}.bam O=${outdir}/align/${prefix}.insert_size_metrics.txt H=${outdir}/align/${prefix}.insert_size_histogram.pdf M=0.5
java -jar ${Picarddir}/picard.jar MarkDuplicates I=${outdir}/align/${prefix}.bam O=${outdir}/align/${prefix}.rmd.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${outdir}/align/${prefix}.metrics REMOVE_DUPLICATES=false && samtools flagstat ${outdir}/align/${prefix}.rmd.bam > ${outdir}/align/${prefix}.rmd.bam.stat

python ${Tooldir}/UniqFileBam.py -t bwa-mem -i ${outdir}/align/${prefix}.rmd.bam -o ${outdir}/align/${prefix}.flt.bam -q 20 && samtools index ${outdir}/align/${prefix}.flt.bam && samtools flagstat ${outdir}/align/${prefix}.flt.bam > ${outdir}/align/${prefix}.flt.bam.stat 

# samtools view -o ${outdir}/align/${prefix}.flt.rm.bam ${outdir}/align/${prefix}.flt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 && samtools index ${outdir}/align/${prefix}.flt.rm.bam && samtools flagstat ${outdir}/align/${prefix}.flt.rm.bam > ${outdir}/align/${prefix}.flt.rm.bam.stat

Rscript ${ppddir}/run_spp.R  -c=${outdir}/align/${prefix}.flt.bam -savp=${outdir}/align/${prefix}.flt.bam.spp.pdf -out=${outdir}/align/${prefix}.flt.bam.spp.txt > ${outdir}/align/${prefix}.flt.bam.richtest.log

bamCoverage -o ${outdir}/align/${prefix}.RPKM.bw -bs 5 -b ${outdir}/align/${prefix}.flt.bam -p ${thread} --minMappingQuality 20 --normalizeUsing RPKM && bamCoverage -o ${outdir}/align/${prefix}.bw -bs 5 -b ${outdir}/align/${prefix}.flt.bam -p ${thread} --minMappingQuality 20

mkdir -p ${outdir}/peak
macs2 callpeak -t ${outdir}/align/${prefix}.flt.bam -c ${input.bam} -f BAM -g ${GENOME_LEN} -n ${outdir}/peak/${prefix}.BAM.narrow --trackline -B --verbose 3 --SPMR --keep-dup all && macs2 callpeak -t ${outdir}/align/${prefix}.flt.bam -c ${input.bam} -f BAM -g ${GENOME_LEN} -n ${outdir}/peak/${prefix}.BAM.broad --trackline -B --verbose 3 --SPMR --broad --keep-dup all
macs2 callpeak -t ${outdir}/align/${prefix}.flt.bam -c ${input.bam} -f BAMPE -g ${GENOME_LEN} -n ${outdir}/peak/${prefix}.BAMPE.narrow --trackline -B --verbose 3 --SPMR --keep-dup all && macs2 callpeak -t ${outdir}/align/${prefix}.flt.bam -c ${input.bam} -f BAMPE -g ${GENOME_LEN} -n ${outdir}/peak/${prefix}.BAMPE.broad --trackline -B --verbose 3 --SPMR --broad --keep-dup all

bedtools multicov -bams ${outdir}/align/${prefix}.flt.bam -bed ${outdir}/peak/${prefix}.BAM.narrow_peaks.narrowPeak > ${outdir}/peak/${prefix}.BAM.narrow_peaks.narrowPeak.count
bedtools multicov -bams ${outdir}/align/${prefix}.flt.bam -bed ${outdir}/peak/${prefix}.BAM.broad_peaks.broadPeak > ${outdir}/peak/${prefix}.BAM.broad_peaks.broadPeak.count
bedtools multicov -bams ${outdir}/align/${prefix}.flt.bam -bed ${outdir}/peak/${prefix}.BAMPE.narrow_peaks.narrowPeak > ${outdir}/peak/${prefix}.BAMPE.narrow_peaks.narrowPeak.count
bedtools multicov -bams ${outdir}/align/${prefix}.flt.bam -bed ${outdir}/peak/${prefix}.BAMPE.broad_peaks.broadPeak > ${outdir}/peak/${prefix}.BAMPE.broad_peaks.broadPeak.count

# multiqc -o ${outdir}/all_combine_raw_multi ${outdir}/rawqc
# multiqc -o ${outdir}/all_combine_filter.multi ${outdir}/trim




#### Trimmomatic
## Trimmomatic是专门用来过滤数据质量的软件，在建库是，根据不同的建库手段，不同数据所携带的接头序列是不一样的，这里给出一般的情况
RNA-Seq、ChIP-Seq、WGS、Hi-C：Illumina Universal Adapter（adapters/TruSeq3-PE-2.fa:2:30:7:8:true or adapters/TruSeq3-PE-2.fa:2:30:10:8:true）
ATAC-Seq、ChIA-PET、ChRD-PET：Nextera Transposase Sequence

单双端使用：
SE ${IN1} ${outdir}/trim/${prefix}_R1.fq.gz
PE ${IN1} ${IN2} ${outdir}/trim/${prefix}_R1.fq.gz ${outdir}/trim/${prefix}_R1.unpair.fq.gz ${outdir}/trim/${prefix}_R2.fq.gz ${outdir}/trim/${prefix}_R2.unpair.fq.gz

#### bwa
## bwa专门用来比对DNA来源的二代测序数据，分成4种情况
1. 长双端数据（一般指测序reads大于70bp）
bwa mem -t ${thread} $index ${outdir}/trim/${i}_R1.fq.gz ${outdir}/trim/${i}_R2.fq.gz -R "@RG\tID:${prefix}\tSM:${prefix}\tLB:${species}\tPL:illumina\tPU:run" | samtools view -Sb - | samtools sort -@ ${thread} -o ${outdir}/align/${prefix}.bam -
2. 短双端数据（一般指测序reads小于50bp）
3. 长单端数据（一般指测序reads大于70bp）
4. 短单端数据（一般指测序reads小于50bp）
bwa aln -t ${thread} $index -f ${outdir}/align/${prefix}.sai ${outdir}/trim/${i}_R2.fq.gz && bwa samse -r "@RG\tID:${prefix}\tSM:${prefix}\tLB:${species}\tPL:illumina\tPU:run" $index ${outdir}/align/${prefix}.sai ${outdir}/trim/${i}_R2.fq.gz | samtools view -Sb - | samtools sort -@ {thread} -o ${outdir}/align/${prefix}.bam -

#### idr
## idr专门用来合并两个ChIP-Seq的peak结果。
idr --samples narrowPeak1 narrowPeak2 --input-file-type narrowPeak --output-file idrValues.peak --output-file-type narrowPeak --plot idrValues.peak.png --verbose