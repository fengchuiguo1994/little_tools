```
#### prepare data
python getRead1M.py 10xARC/outs/gex_possorted_bam.bam ../SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq
python getRead1M.py 10xARC/outs/gex_possorted_bam.bam ../SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq
python getRead1M.py 10xARC/outs/atac_possorted_bam.bam ../SCG0074-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R1_001.fastq
python getRead1M.py 10xARC/outs/atac_possorted_bam.bam ../SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq
python getRead1M.py 10xARC/outs/atac_possorted_bam.bam ../SCG0074-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R3_001.fastq
gzip SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq
gzip SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq
gzip SCG0074-15872_SI-NA-D6_S5_L003_R1_001.fastq
gzip SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq
gzip SCG0074-15872_SI-NA-D6_S5_L003_R3_001.fastq

#### cellranger-arc process
cellranger-arc count --reference genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --id 10xARC --libraries=SCG0074.arc.csv --localcores=10 --jobmode=local

#### STARsolo
singularity exec ~/Tools/STAR.sif STAR --genomeDir genome/stargenome/hg38_index --outFileNamePrefix SCG0081_STARsoloCellRanger4 --readFilesCommand zcat --readFilesIn SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 737K-arc-v1.txt --soloCellFilter EmptyDrops_CR --runThreadN 10 --outSAMattributes CB UB --outSAMtype BAM SortedByCoordinate --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR

#### chromap
paste 737K-arc-v1.txt ATAC.737K-arc-v1.txt > 737K-arc_RNA_ATAC_barcodeTranslation.tsv
chromap -t 10 --preset atac -x /data/home/ruanlab/huangxingyu/Tools/cellranger-arc-2.0.2/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/chromap/genonme.index -r /data/home/ruanlab/huangxingyu/Tools/cellranger-arc-2.0.2/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa -1 SCG0074-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz -2 SCG0074-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz -b SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz --barcode-whitelist ATAC.737K-arc-v1.txt --SAM -o chromapSCG0074_fragments.ts.sam --read-format bc:8:23:- --barcode-translate 737K-arc_RNA_ATAC_barcodeTranslation.tsv >chromap.log 2>chromap.err
```