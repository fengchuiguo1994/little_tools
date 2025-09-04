- ChIAPET.sh 
分析ChIAPET/ChIATAC数据，识别连接子，比对，数据拆分为UU、Ux、none，去除F2048，去重，识别峰，识别环

有关ChIATAC分析核心软件CPU的一些说明 (ChIATAC-CPU说明.pptx)

#### chiapet2hic
```
# chiapet_tools流程的
i=SCG0192
awk -v OFS="\t" '{print 0,$1,$2,0,0,$4,$5,1}' ${i}.bedpe.selected.pet.txt | sort --parallel=5 -k2,2 -k6,6 -k3,3n -k7,7n > ${i}.chiapet.ValidPairs
java -jar /data/home/ruanlab/huangxingyu/20230615/xiaoqin/230602C-S-YXH-Z01/juicer_tools.1.8.9_jcuda.0.8.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000 ${i}.chiapet.ValidPairs ${i}.chiapet.ValidPairs.hic mm10

# 基于cpu的流程
/data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/pairix_src/util/bam2pairs/bam2pairs -c /data/home/ruanlab/huangxingyu/Genome/hg38/hg38.chrom.sizes /data/home/ruanlab/huangxingyu/dataTest/ChIA-PETTEST/result/step2/LHH0135.singlelinker.paired.UU.nr.bam LHH0135
java -Xmx30g -jar /data/home/ruanlab/huangxingyu/Tools/ChIA-PIPE-master/util/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 LHH0135.bsorted.pairs.gz ChIA-PET_hg38_HG00513_CTCF_LHH0135_miseq_pairs.hic /data/home/ruanlab/huangxingyu/Genome/hg38/hg38.chrom.sizes
```

#### chiapet数据过滤
```
java -jar /data/home/ruanlab/huangxingyu/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 -phred33 LHH0135_GT20-15524_GGACTCCT-ATAGAGAG_S1_L002_R1_001.fastq.gz LHH0135_GT20-15524_GGACTCCT-ATAGAGAG_S1_L002_R2_001.fastq.gz LHH0135_R1.fq.gz LHH0135_R1.un.fq.gz LHH0135_R2.fq.gz LHH0135_R2.un.fq.gz LEADING:10 TRAILING:10 ILLUMINACLIP:/data/home/ruanlab/huangxingyu/Tools/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:7:8:true SLIDINGWINDOW:4:15 MINLEN:50
```

#### chiapet连接两端的coverage
```
awk '$2<$5 {print $1 "\t" $2 "\t" $5}; $2>$5 {print $1 "\t" $5 "\t" $2}' RMCP-101.ipet | sort-k 1.4n -k2n -k3n | awk '{print $0 "\t" "IPET"NR}' | bedToBam -i stdin -g chrom.sizes > RMCP-101.ipet.sorted.bam
samtools index RMCP-101.ipet.sorted.bam
multiBamSummary bins --bamfiles *.ipet.sorted.bam --minMappingQuality 30 -out RMCP-101 RMCP-093 ipet.npz --outRawCounts RMCP-101 RMCP-093ipet.tab
