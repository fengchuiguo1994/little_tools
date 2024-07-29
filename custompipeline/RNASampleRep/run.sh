# cp /data/home/ruanlab/hanjin/RNA-seq/HCRS000?/2_map/2_clean/HCRS000?.clean.bam* .

for i in HCRS0001 HCRS0002 HCRS0003
do
    # bamCoverage -b ${i}.clean.bam -o ${i}.bw -p 10
    # bamCoverage -b ${i}.clean.bam -o ${i}.bedgraph -of bedgraph -p 10
    # stringtie ${i}.clean.bam -G ~/Tools/cellranger-7.1.0/genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf -o ${i}.gtf -p 8 -e -b ${i} -A ${i}.abundance
    # awk '{if(NR==FNR){aa[$3"\t"$5"\t"$6"\t"$1]=$0}else{if($1"\t"$2+1"\t"$3"\t"$4 in aa){print aa[$1"\t"$2+1"\t"$3"\t"$4];}else{print $4"\t0\t0\t0\t0\t0\t0\t0\t0"}}}' ${i}.abundance /data/home/ruanlab/huangxingyu/Tools/cellranger-7.1.0/genome/refdata-gex-GRCh38-2020-A/genes/hg38.info > ${i}.flt.abundance
    sleep 0.01
done

# multiBigwigSummary bins -b HCRS0001.bw HCRS0002.bw HCRS0003.bw -out HCRS.results.npz --labels HCRS0001 HCRS0002 HCRS0003 --numberOfProcessors 10 --outRawCounts HCRS.results.RawCounts
# grep -v -e chrUn_ -e _random -e _alt HCRS.results.RawCounts > HCRS.results.flt.RawCounts

# featureCounts -a ~/Tools/cellranger-7.1.0/genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf -o HCRS.gene.count -t exon -g gene_id -p -T 1 HCRS0001.clean.bam HCRS0002.clean.bam HCRS0003.clean.bam



# paste <(cut -f 1,9 HCRS0001.flt.abundance) <(cut -f 9 HCRS0002.flt.abundance) <(cut -f 9 HCRS0003.flt.abundance) > HCRS.tpm
# paste <(cut -f 1,8 HCRS0001.flt.abundance) <(cut -f 8 HCRS0002.flt.abundance) <(cut -f 8 HCRS0003.flt.abundance) > HCRS.fpkm


