# 文献检索
[scihub](https://www.scihub.net.cn/)<br/>
[scihub2](https://tool.yovisun.com/scihub/)<br/>
[wiki百科](https://zh.wikipedia.org/wiki/Wikipedia:%E9%A6%96%E9%A1%B5)<br/>

# RNA 数据库条目
### snoRNA/snRNA database
刚开始的时候，snoRNA和snRNA没有明显的差别，可以看到很多早期文献还在称U3为snRNA。所以合并了snRNA/snoRNA。<br/>
[酵母的snoRNA数据库](https://people.biochem.umass.edu/fournierlab/snornadb/mastertable.php)<br/>
[人的snoRNA数据库](https://www-snorna.biotoul.fr/getseq.php)。可以从这个网站下载snoRNA的序列。<br/>

[植物单细胞RNA-seq数据库](http://ibi.zju.edu.cn/plantscrnadb/index.php)<br/>

### ncRNA(有冗余分类)
从Rfam中提取ncRNA的序列。
```
mkdir ~/rfam && cd ~/rfam
curl -C -O http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF0[0001-4222].fa.gz
gunzip *.gz
```
[植物circRNA](http://ibi.zju.edu.cn/index.html/sever.html)<br/>
[GREAT](http://great.stanford.edu/public/html/)，根据附近的基因的注释，对非编码基因区域进行生物学意义注释。<br/>
[AnnoLnc2](http://annolnc.gao-lab.org/index.php)。人和小鼠lncRNA的注释，包括lncRNA预测，相关miRNA，保守性，二级结构，GO注释等等。[使用 RepeatMasker 来预测的序列当中是否有重复序列; 通过ViennaRNA (http://rna.tbi.univie.ac.at/) 数据库来预测lncRNA的二级结构；数据库使用了GETx数据库里面的正常组织、CCLE里面的癌症细胞系以及ENCODE数据库里面的数据来进行查看基因的表达情况；通过比较核/胞质表达来确定这个lncRNA主要是在哪个地方表达；使用GTRD来预测lncRNA的可能收到的转录因子调控作用，同时使用TargetScan来预测其miRNA调控的作用；通过GWAS数据库来寻找影响这个lncRNA的SNP，进一步的通过eQTL来评价哪些SNP对于这个lncRNA的表达有影响，这个分析的主要数据来自于GETx；由于使用的RNA-seq的数据，所以就可以看lncRNA的表达和哪些基因存在共表达关系；使用了目前发表的GEO上面的CLIP-seq的数据来进行分析，对于GEO里面没有的蛋白数据，数据库使用lncPro数据库来进行预测。所以在结果当中就包括两个部分，一个是lncPro数据库的结果，另外一个则是CLIP-seq分析的结果；预测这个lncRNA的功能了。由于lncRNA本身是不会编码蛋白来发挥作用的，所以主要是通过其相互作用的基因来预测这个lncRNA的功能，这个数据库主要预测了lncRNA本身GO分析的功能；通过phyloFit来比较物种之间的进化关系](https://www.sci666.com.cn/66870.html)<br/>
[人miRNA的靶基因](http://www.targetscan.org/vert_80/)<br/>
[实验室里所有的数据库，包括plantTFDB，PlantRegMap，CPC，AnnoLnc2](http://plantregmap.cbi.pku.edu.cn/)<br/>

### miRNA
[TargetScan](https://www.targetscan.org/vert_72/)是一款预测miRNA结合位点的软件，对于哺乳动物中miRNA结合位点预测的效果非常好。在预测miRNA靶基因之前，首先需要确定转录本的3’UTR区域，TargetScan数据库通过一种名为3P-seq的测序技术，确定转录本对应的3’UTR区（哺乳动物中的miRNA通过结合转录本序列的3’UTR区，从而发挥转录后调控作用），并且结合该技术的分析结果和NCBI中已有的3’UTR注释，提供一个综合的3’UTR区序列。[来源于知乎](https://zhuanlan.zhihu.com/p/325673305)<br/>
[miRanda](http://cbio.mskcc.org/miRNA2003/miranda.html)，[文档](https://bioinformaticsworkbook.org/dataAnalysis/SmallRNA/Miranda_miRNA_Target_Prediction.html#gsc.tab=0)<br/>


# RNA 二级结构
[RNArchitecture](http://genesilico.pl/RNArchitecture/family/MIR807/secondarystructure)<br/>
[ViennaRNA预测RNA的二级结构](https://www.tbi.univie.ac.at/RNA/)<br/>
[预测RNA的二级结构](http://rna.tbi.univie.ac.at/)<br/>

# repeat 数据库条目（有些repeat RNA也属于该条目）
注释条目说明
```
Short interspersed nuclear elements (SINE), which include ALUs
Long interspersed nuclear elements (LINE)
Long terminal repeat elements (LTR), which include retroposons
DNA repeat elements (DNA)
Simple repeats (micro-satellites)
Low complexity repeats
Satellite repeats
RNA repeats (including RNA, tRNA, rRNA, snRNA, scRNA, srpRNA)
Other repeats, which includes class RC (Rolling Circle)
Unknown
```
### repeat database
[RepBase](https://www.girinst.org/server/RepBase/)是一个整理的repeat数据库，但是使用需要支付昂贵的费用。在repeatmasker网站上公开了两个版本的[RepBase](http://repeatmasker.org/libraries/)：RepBaseRepeatMaskerEdition-20181026.tar.gz和RepeatMaskerMetaData-20170127.tar.gz。<br/>
[Dfam](https://www.dfam.org/home), [Transposable Element DNA sequence alignments, hidden Markov Models (HMMs)](https://www.dfam.org/releases/Dfam_3.5/annotations/hg38/)<br/>
### repeatmask注释条目
[RMGenomicDatasets](http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html)年代久远而且和UCSC下载的对不上。UCSC的用的是20130422版。<br/>
[repeatmask的说明](http://www.repeatmasker.org/faq.html).<br/>

# genome & annotation 数据库条目
### 水稻基因组
[ZS97 & MH63](https://rice.hzau.edu.cn/cgi-bin/rice_rs3/download_ext)，缺点是不含线粒体和叶绿体<br/>
[日本晴 Oryza](https://rapdb.dna.affrc.go.jp/download/irgsp1.html)，含有线粒体，叶绿体<br/>
[UCSC基因组来源](http://genome.ucsc.edu/goldenPath/credits.html)<br/>
[USCS收录的人类基因组注释](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/)，以及一个[更新版](http://hgdownload.soe.ucsc.edu/goldenPath/archive/hg38/)<br/>
[human基因注释和功能注释 genecards](https://www.genecards.org/)<br/>
[human基因注释和功能注释 HGNC](https://www.genenames.org/)<br/>
[植物数据库导航](https://mp.weixin.qq.com/s?__biz=MzA5OTI1MzM4Mw==&mid=2247485282&idx=1&sn=8c14fa485455a32bea5e2a269ceedcb8&chksm=90846a8aa7f3e39cf9a23dc4c8ea62ef696a1dce6f2d4648d316b7c3ed0c7aef55c788f17006&mpshare=1&scene=1&srcid=1217ogKVlbYaJWeV5ukQtEoH&sharer_sharetime=1608189858165&sharer_shareid=e80cfb57b06d9707ef6ec23eb958dc90&key=48ac2953e3e884862a8ced436c4c1dbb6998cdb350c32cd5451da0a0998d5a2ee85b1e296bf8b337fa422b1c689694c7177a7e02c63a8bb27e145b7f38a057f19b17c98eb50ce4ba5877c19358ea0b7efb9d7d86c679562faf47652a62bfcb90f81d20c1a0f34e7dea4a38da8557a715b3318c527677ccf0a6005d5b80886be2&ascene=1&uin=MjE5OTExMDU3Ng%3D%3D&devicetype=Windows+10+x64&version=6300002f&lang=zh_CN&exportkey=AUzUCLzZ8FPBQMQV%2Bz3T1Dw%3D&pass_ticket=epI%2Bxd1QO5YvutZSIErkHCb9QjcGayq3GbcmlXIAGa2ek8%2BUZn8qrIlkHP6fLfjW&wx_header=0)<br/>

### 玉米基因组
[B73](https://www.maizegdb.org/assembly#stockInfo)<br/>

### 小鼠基因组
[ChromHMM注释](https://github.com/gireeshkbogu/chromatin_states_chromHMM_mm9)<br/>

### 人类基因组
[gwas_catalog]() gwas数据库<br/>
[GRC](https://www.ncbi.nlm.nih.gov/grc)<br/>
[ChromHMM注释](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html)结果<br/>
[ChromHMM解释](https://pubs.broadinstitute.org/mammals/haploreg/documentation_v2.html)<br/>
```
STATE NO.	MNEMONIC	DESCRIPTION	COLOR NAME	COLOR CODE
1	TssA	Active TSS	Red	255,0,0
2	TssAFlnk	Flanking Active TSS	Orange Red	255,69,0
3	TxFlnk	Transcr. at gene 5' and 3'	LimeGreen	50,205,50
4	Tx	Strong transcription	Green	0,128,0
5	TxWk	Weak transcription	DarkGreen	0,100,0
6	EnhG	Genic enhancers	GreenYellow	194,225,5
7	Enh	Enhancers	Yellow	255,255,0
8	ZNF/Rpts	ZNF genes & repeats	Medium Aquamarine	102,205,170
9	Het	Heterochromatin	PaleTurquoise	138,145,208
10	TssBiv	Bivalent/Poised TSS	IndianRed	205,92,92
11	BivFlnk	Flanking Bivalent TSS/Enh	DarkSalmon	233,150,122
12	EnhBiv	Bivalent Enhancer	DarkKhaki	189,183,107
13	ReprPC	Repressed PolyComb	Silver	128,128,128
14	ReprPCWk	Weak Repressed PolyComb	Gainsboro	192,192,192
15	Quies	Quiescent/Low	White	255,255,255

STATE NO.	MNEMONIC	DESCRIPTION	COLOR NAME	COLOR CODE
1	TssA	Active TSS	Red	255,0,0
2	TssFlnk	Flanking TSS	Orange Red	255,69,0
3	TssFlnkU	Flanking TSS Upstream	Orange Red	255,69,0
4	TssFlnkD	Flanking TSS Downstream	Orange Red	255,69,0
5	Tx	Strong transcription	Green	0,128,0
6	TxWk	Weak transcription	DarkGreen	0,100,0
7	EnhG1	Genic enhancer1	GreenYellow	194,225,5
8	EnhG2	Genic enhancer2	GreenYellow	194,225,5
9	EnhA1	Active Enhancer 1	Orange	255,195,77
10	EnhA2	Active Enhancer 2	Orange	255,195,77
11	EnhWk	Weak Enhancer	Yellow	255,255,0
12	ZNF/Rpts	ZNF genes & repeats	Medium Aquamarine	102,205,170
13	Het	Heterochromatin	PaleTurquoise	138,145,208
14	TssBiv	Bivalent/Poised TSS	IndianRed	205,92,92
15	EnhBiv	Bivalent Enhancer	DarkKhaki	189,183,107
16	ReprPC	Repressed PolyComb	Silver	128,128,128
17	ReprPCWk	Weak Repressed PolyComb	Gainsboro	192,192,192
18	Quies	Quiescent/Low	White	255,255,255
```
[circRNADb](http://reprod.njmu.edu.cn/cgi-bin/circrnadb/circRNADb.php)环状RNA<br/>
```
cd ~/reference
mkdir -p genome/hg19  && cd genome/hg19 
nohup wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz &
tar zvfx chromFa.tar.gz
cat *.fa > hg19.fa
rm chr*.fa
 
cd ~/reference
mkdir -p genome/hg38  && cd genome/hg38 
nohup wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz  &
 
cd ~/reference
mkdir -p  genome/mm10  && cd genome/mm10 
nohup wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz  &
tar zvfx chromFa.tar.gz
cat *.fa > mm10.fa
rm chr*.fa
 
cd ~/biosoft/RNA-SeQC
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/ThousandReads.bam
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/gencode.v7.annotation_goodContig.gtf.gz
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/Homo_sapiens_assembly19.fasta.gz
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/Homo_sapiens_assembly19.other.tar.gz
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/gencode.v7.gc.txt
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/rRNA.tar.gz
 
cd ~/reference
mkdir -p index/bowtie && cd index/bowtie 
nohup time ~/biosoft/bowtie/bowtie2-2.2.9/bowtie2-build  ~/reference/genome/hg19/hg19.fa  ~/reference/index/bowtie/hg19 1>hg19.bowtie_index.log 2>&1 &
nohup time ~/biosoft/bowtie/bowtie2-2.2.9/bowtie2-build  ~/reference/genome/hg38/hg38.fa  ~/reference/index/bowtie/hg38 1>hg38.bowtie_index.log 2>&1 &
nohup time ~/biosoft/bowtie/bowtie2-2.2.9/bowtie2-build  ~/reference/genome/mm10/mm10.fa  ~/reference/index/bowtie/mm10 1>mm10.bowtie_index.log 2>&1 &
  
cd ~/reference
mkdir -p index/bwa && cd index/bwa 
nohup time ~/biosoft/bwa/bwa-0.7.15/bwa index   -a bwtsw   -p ~/reference/index/bwa/hg19  ~/reference/genome/hg19/hg19.fa 1>hg19.bwa_index.log 2>&1   &
nohup time ~/biosoft/bwa/bwa-0.7.15/bwa index   -a bwtsw   -p ~/reference/index/bwa/hg38  ~/reference/genome/hg38/hg38.fa 1>hg38.bwa_index.log 2>&1   &
nohup time ~/biosoft/bwa/bwa-0.7.15/bwa index   -a bwtsw   -p ~/reference/index/bwa/mm10  ~/reference/genome/mm10/mm10.fa 1>mm10.bwa_index.log 2>&1   &
  
cd ~/reference
mkdir -p index/hisat && cd index/hisat 
nohup wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg19.tar.gz  &
nohup wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz  &
nohup wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38.tar.gz &
nohup wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz  &
tar zxvf hg19.tar.gz
tar zxvf grcm38.tar.gz
tar zxvf hg38.tar.gz
tar zxvf mm10.tar.gz 
  
mkdir -p ~/annotation/variation/human/ExAC
cd ~/annotation/variation/human/ExAC
## http://exac.broadinstitute.org/
## ftp://ftp.broadinstitute.org/pub/ExAC_release/current
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/current/ExAC.r0.3.1.sites.vep.vcf.gz.tbi 
nohup wget ftp://ftp.broadinstitute.org/pub/ExAC_release/current/ExAC.r0.3.1.sites.vep.vcf.gz &
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/current/cnv/exac-final-cnv.gene.scores071316 
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/current/cnv/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed
 
mkdir -p ~/annotation/variation/human/dbSNP
cd ~/annotation/variation/human/dbSNP
## https://www.ncbi.nlm.nih.gov/projects/SNP/
## ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh38p2/
## ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/
nohup wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz &
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz.tbi 
 
mkdir -p ~/annotation/variation/human/1000genomes
cd ~/annotation/variation/human/1000genomes 
## ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ 
nohup wget  -c -r -nd -np -k -L -p  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502 &
 
mkdir -p ~/annotation/variation/human/cosmic
cd ~/annotation/variation/human/cosmic
## we need to register before we can download this file. 
 
mkdir -p ~/annotation/variation/human/ESP6500
cd ~/annotation/variation/human/ESP6500
# http://evs.gs.washington.edu/EVS/
nohup wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz & 
 
mkdir -p ~/annotation/variation/human/UK10K
cd ~/annotation/variation/human/UK10K
# http://www.uk10k.org/
nohup wget ftp://ngs.sanger.ac.uk/production/uk10k/UK10K_COHORT/REL-2012-06-02/UK10K_COHORT.20160215.sites.vcf.gz & 
 
mkdir -p ~/annotation/variation/human/gonl
cd ~/annotation/variation/human/gonl
## http://www.nlgenome.nl/search/
## https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5/
nohup wget  -c -r -nd -np -k -L -p  https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5  &
 
mkdir -p ~/annotation/variation/human/omin
cd ~/annotation/variation/human/omin
 
mkdir -p ~/annotation/variation/human/GWAS
cd ~/annotation/variation/human/GWAS
 
mkdir -p ~/annotation/variation/human/hapmap
cd ~/annotation/variation/human/hapmap
# ftp://ftp.ncbi.nlm.nih.gov/hapmap/
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/relationships_w_pops_051208.txt 
nohup wget -c -r -np -k -L -p  -nd -A.gz ftp://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/hapmap3_reformatted &
# ftp://ftp.hgsc.bcm.tmc.edu/pub/data/HapMap3-ENCODE/ENCODE3/ENCODE3v1/
wget ftp://ftp.hgsc.bcm.tmc.edu/pub/data/HapMap3-ENCODE/ENCODE3/ENCODE3v1/bcm-encode3-QC.txt 
wget ftp://ftp.hgsc.bcm.tmc.edu/pub/data/HapMap3-ENCODE/ENCODE3/ENCODE3v1/bcm-encode3-submission.txt.gz
 
## 1 million single nucleotide polymorphisms (SNPs) for DNA samples from each of the three ethnic groups in Singapore – Chinese, Malays and Indians.
## The Affymetrix Genome-Wide Human SNP Array 6.0   && The Illumina Human1M single BeadChip 
## http://www.statgen.nus.edu.sg/~SGVP/
## http://www.statgen.nus.edu.sg/~SGVP/singhap/files-website/samples-information.txt
# http://www.statgen.nus.edu.sg/~SGVP/singhap/files-website/genotypes/2009-01-30/QC/
 
## Singapore Sequencing Malay Project (SSMP) 
mkdir -p ~/annotation/variation/human/SSMP
cd ~/annotation/variation/human/SSMP
## http://www.statgen.nus.edu.sg/~SSMP/
## http://www.statgen.nus.edu.sg/~SSMP/download/vcf/2012_05 
 
## Singapore Sequencing Indian Project (SSIP) 
mkdir -p ~/annotation/variation/human/SSIP
cd ~/annotation/variation/human/SSIP
# http://www.statgen.nus.edu.sg/~SSIP/
## http://www.statgen.nus.edu.sg/~SSIP/download/vcf/dataFreeze_Feb2013
 
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz 
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.chr.gtf.gz 
 
mkdir -p ~/reference/gtf/gencode
cd  ~/reference/gtf/gencode
## https://www.gencodegenes.org/releases/current.html
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.2wayconspseudos.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.long_noncoding_RNAs.gtf.gz 
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.polyAs.gtf.gz 
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz 
## https://www.gencodegenes.org/releases/25lift37.html 
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz 
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.HGNC.gz 
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.EntrezGene.gz 
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.metadata.RefSeq.gz 
 
mkdir -p ~/reference/gtf/ensembl/homo_sapiens_86
cd  ~/reference/gtf/ensembl/homo_sapiens_86
## http://asia.ensembl.org/info/data/ftp/index.html
 
cd ~/reference
mkdir -p  genome/human_g1k_v37  && cd genome/human_g1k_v37
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/ 
nohup wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz  &
gunzip human_g1k_v37.fasta.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/README.human_g1k_v37.fasta.txt
java -jar ~/biosoft/picardtools/picard-tools-1.119/CreateSequenceDictionary.jar R=human_g1k_v37.fasta O=human_g1k_v37.dict
 
## ftp://ftp.broadinstitute.org/bundle/b37/
mkdir -p ~/annotation/GATK
cd ~/annotation/variation/GATK
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz 
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz 
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz 
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz 
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.idx.gz
gunzip 1000G_phase1.indels.b37.vcf.idx.gz
gunzip 1000G_phase1.indels.b37.vcf.gz
  
mkdir -p  ~/institute/ENSEMBL/gtf
cd  ~/institute/ENSEMBL/gtf
wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.chr.gtf.gz 
wget ftp://ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.chr.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-87/gtf/danio_rerio/Danio_rerio.GRCz10.87.chr.gtf.gz

cd ~/institute/TCGA/firehose
## https://gdac.broadinstitute.org/
wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ACC/20160128/gdac.broadinstitute.org_ACC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0.tar.gz  -O ACC.gistic.seg.tar.gz
wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ACC/20160128/gdac.broadinstitute.org_ACC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2016012800.0.0.tar.gz  -O ACC.raw.seg.tar.gz 
wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ACC/20160128/gdac.broadinstitute.org_ACC.Mutation_Packager_Calls.Level_3.2016012800.0.0.tar.gz -O ACC.maf.tar.gz
wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ACC/20160128/gdac.broadinstitute.org_ACC.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0.tar.gz -O ACC.maf.anno.tar.gz
```

### 植物基因组
[PlantGDB 植物基因组整合](http://www.plantgdb.org/)

# 常用在线小工具
### 比对
[多序列比对](https://www.novopro.cn/tools/muscle.html)<br/>
[DeepL翻译器](https://www.deepl.com/translator)<br/>
[科研者之家](https://www.home-for-researchers.com/static/index.html#/)<br/>
网络图绘制，可以用[cytoscape](https://cytoscape.org/)和[Gephi](https://gephi.org/)<br/>
InterProscan 输出格式
```
1. 蛋白质接入号	Protein Accession (e.g. P51587)
2. 序列的 MD5 值	Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
3. 序列长度	Sequence Length (e.g. 3418)
4. 不同分析方案	Analysis (e.g. Pfam / PRINTS / Gene3D)
5. 签名号	Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
6. 签名描述	Signature Description (e.g. BRCA2 repeat profile)
7. 起始位置	Start location
8. 终止位置	Stop location
9. 得分	Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
10. 状态	Status - is the status of the match (T: true)
11. 运行日期	Date - is the date of the run
12... 其他 	(InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
	(InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
	(GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
	(Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
```
[InterProscan的interpro2go注释文件](ftp.ebi.ac.uk/pub/databases/interpro/interpro2go)，[所有数据](http://ftp.ebi.ac.uk/pub/databases/interpro/)
# 林奈分类以及一些常见简写

简写 | 全称 | 翻译 | 发现的场景
----|----|----|----
CB  | Caenorhabditis briggsae | 广杆属线虫 | repbase中查阅
ECa | Eucalyptus camaldulensis| 赤桉; 桉树林多为赤桉 | repbase中查阅
EC  | Equus caballus | 马 | repbase中查阅（少量）
ML  | Myotis lucifugus|蝙蝠|repbase中查阅
CR  | Chlamydomonas reinhardtii | 莱茵衣藻 | repbase中查阅
SP  | Strongylocentrotus | 海胆 | repbase中查阅

# 植物数据库
### GWAS数据库
[gwas Atlas](https://ngdc.cncb.ac.cn/gwas/)
### 水稻
[funRiceGenes 水稻基因功能基因名注释](https://funricegenes.github.io/)<br/>
[rap-db 水稻基因组](https://rapdb.dna.affrc.go.jp/index.html)<br/>
[Oryzabase 水稻注释](https://shigen.nig.ac.jp/rice/oryzabase/)<br/>
[RIGW ZS97和MH63的详细注释](http://rice.hzau.edu.cn/rice_rs3/)<br/>
[Gene lncRNA circRNA数据库](http://ic4r.org/)<br/>
[综合数据库字典](https://www.hsls.pitt.edu/obrc/index.php?page=rice)<br/>

### 复合植物
[大豆 拟南芥 水稻 玉米RNA-seq数据库](http://ipf.sustech.edu.cn/pub/athrdb/)<br/>
[植物repeat 保守性 开放区 GO注释 转录因子等，个人只用过保守性的结果](http://plantregmap.gao-lab.org/download.php#alignment-conservation)<br/>
[PlantTFDB 植物转录因子数据库，跟上面一起的](http://planttfdb.gao-lab.org/)<br/>
[植物复合数据库](https://mp.weixin.qq.com/s?__biz=MzA5OTI1MzM4Mw==&mid=2247485282&idx=1&sn=8c14fa485455a32bea5e2a269ceedcb8&chksm=90846a8aa7f3e39cf9a23dc4c8ea62ef696a1dce6f2d4648d316b7c3ed0c7aef55c788f17006&mpshare=1&scene=23&srcid=12175Mndpizc8SCj0DaIwUWR&sharer_sharetime=1608189227953&sharer_shareid=fb0a8cc66caa940860c8750158cdb4b4#rd)<br/>
[植物GWAS](https://easygwas.ethz.ch/)<br/>
[动物植物RNA互作数据库](http://rnainter.org/batch_search/)<br/>

# 动物数据库
[人和鼠等多组织表达量数据查询](https://www.ebi.ac.uk/gxa/home#)<br/>
[人的blacklist位置注释](https://www.encodeproject.org/data-standards/reference-sequences/)<br/>
[Cell Marker 细胞标志物](http://xteam.xbio.top/CellMarker/)<br/>
[expasy 查询细胞系详细信息](https://web.expasy.org/cellosaurus/CVCL_0168)<br/>
[RNA的二级结构](https://biocyc.org/gene?orgid=ECOLI&id=gltW-tRNA)<br/>
[人snoRNA数据库 snoRNABase](https://www-snorna.biotoul.fr/)<br/>
[人 GWAS Catalog](https://www.ebi.ac.uk/gwas/search?query=rs13422172)<br/>
[人和小鼠的lncRNA的数据库 AnnoLnc2](http://annolnc.gao-lab.org/)<br/>
[RegulomeDB 人类SNP的gwas和其它关联分析的结果](https://www.regulomedb.org/regulome-search/)<br/>
[HaploReg 人类SNP的gwas和其它关联分析的结果](https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php)<br/>
[RISE 人类RNA互作数据库](http://rise.life.tsinghua.edu.cn/)<br/>

# GO工具
[EMBL-EBI quickGO](https://www.ebi.ac.uk/QuickGO/)<br/>
[AmiGO2](http://amigo.geneontology.org/amigo/term/GO:0005515)<br/>
[ShinyGO](http://bioinformatics.sdstate.edu/go/)<br/>
[metascape](https://metascape.org/gp/index.html#/main/step1)<br/>
[GOplot 画图 R包](https://www.biomart.cn/9594/news/2949857.htm), [参考案例](https://www.jianshu.com/p/48ac98098760)<br/>
```

```

# 印记基因
### [印记基因数据库](https://www.geneimprint.com/site/genes-by-species.Homo+sapiens)


# 蛋白数据库
[nextprot](https://www.nextprot.org/)
# 人类变异的数据库
```
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz.tbi
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3_grch38_pop_stratified_af.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/wgs_calling_regions.hg38.interval_list

ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.idx.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.idx.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz
```

# 杂项
### 分析流程
[4DN组织提供了各种数据分析的流程](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline?redirected_from=%2Fhelp%2Fanalysis-and-visualization%2Fhi_c-processing-pipeline)<br/>
[生物信息分析的流程和工具介绍 字典](http://www.bio-soft.net/rna.html)<br/>
[RNA结构预测](https://www.cnblogs.com/jcf666/p/13753415.html)<br/>