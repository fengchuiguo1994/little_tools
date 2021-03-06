# some little tools for bioinformatics

python在服务器上画图，开头添加
```python
import matplotlib
matplotlib.use('agg')
```

### cal-RNA-file.py
Statistical RNA-Seq data alignment results(hisat2/star/magic-blast/tophat2)<br/>
python ~/cal-RNA-file.py input.sam/bam

### flt_bam.py
split hisat2 alignment file to un/uniq/multi-map file
python flt_bam.py input.bam/sam output_prefix

### overlap.pl
Calculate the distance between two intervals(chr1   start1  end1    chr2    start2  end2)<br/>
perl overlap.pl bedpe.bed

### compare_exon_intron_coverage.py
if we don't know what type of file(DNA or RNA),wo can detect it through comparing exon and intron coverage,if it RNA file,the exon have higher coverage <br/>
python compare_exon_intron_coverage.py

### reptile.136book.py
using urllib and bs4 to get the story from 136book.com<br/>
python reptile.136book.py

### word_cloud_test.py
using jieba and pyecharts to generate the word cloud<br/>
python word_cloud_test.py

### count2FPKM.py
using htseq-count output and gene exon length to calculator the [FPKM and FPKM-UQ](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)<br/>
python count2FPKM htseq.count exon_gene_length.txt out.txt<br/>
different with TCGA:<br/>
- the total reads just use protein_coding gene to calculate the total reads and Q75 reads

### annopeak.py
using bedtools to annotation peak<br/>
python annopeak.py peak/allrich-06.DNA.promoter.overlap.bed peak/allrich-06.DNA.genebody.overlap.bed peak/allrich-06.DNA.exon.overlap.bed peak/allrich-06.DNA.terminal.overlap.bed > peak/allrich-06.DNA.anno.result.txt<br>
future word<br/>
- using pybedtools to do this work.

### parsegtf.py
if use stringtie or cufflinks to assemble,the output gtf file don't contain the gene line.this script is used to add it.<br/>
python parsegtf.py MH63.seedling.nove.annotated.gtf > MH63.seedling.nove.cg.gtf

### pieplot1.r
just plot pie plot
![](img/pie.plot.png)

### combine_fastq.py
对于RNA-Seq,如果是两个重复数据(生物学重复,两个 library)。在第一个文件里的某个 reads 它会比对到基因组的某个位置,然后再将两个文件直接合并,这个同样的 reads 有一定的概率比对到其它的位置,可能这些 RNA-Seq 比对工具(tophat/hisat/star,已验证,有数据集在百度云)会根据 insert size 评估比对的结果导致的。所以合并两个RNA-Seq的时候建议 A-B-A-B-A-B 的形式合并而不是A-A-...-A-A-B-B-...-B的合并。<br/>
python combine.py test1.fq test2.fq test.fq<br/>
待做的<br/>
- 没有考虑数据大小的差异，如果A文件时B文件的n倍大,那么应该时 nA-B-nA-B 的形式。

### bam2fastq.py
对按照name排序的bam文件转化成fastq文件 <br/>
$ `samtools sort -n -o file.sorted.bam file.bam` <br/>
$ `python bam2fq.py ORP-113.multimap.bam ORP-113.multimap`<br/>
将bam/sam文件转化成fastq文件，支持单端比对/双端比对。如果输入的是双端数据,双端数据中成对的结果保存在ORP-113.multimap.N.fq中,如果是单端数据则在ORP-113.multimap.fq中,混合数据也同理

### divide_bam.py
python divide_bam.py -i test.bam -t bwamem -f bam<br/>
将bam/sam文件分成uniq,multi,unmap,lowuniq四类，并且统计<br/>

### UniqFileBam.py
python UniqFileBam.py -t bowtie2 -i mappings.bowtie2.sam -o test2.sam -q 20 -f True<br/>
提取bam/sam文件中的unique比对的结果，支持 bwa-aln bwa-mem bowtie2 hisat/hisat2 STAR (后续用过别的软件持续更新)<br/>

### htseq.py
python htseq.py<br/>
修改自htseq-count,可以同时计算sense(yes),no-sense(no),anti-sense(reverse)的结果,目前只适合单端数据,日后也没有添加双端的计划,STAR的count的结果默认就有三种<br/>

### split_bam_SE.py
python split_bam_SE.py input.bam input <br/>
根据指定的文件行数分割文件，暂时没有提供接口，目前都是切分1000000左右<br/>

### split_bam_PE.py
python split_bam_SE.py input.bam input <br/>
类上的split_bam_SE.py，根据指定的文件行数分割文件，暂时没有提供接口，目前都是切分1000000左右<br/>

### add_count.py
python add_count.py input 1 39 > total.count <br/>
通过上面的 split_bam 将文件分割，然后用 htseq-count 统计基因的表达数量，通过该软件合并起来<br/>

### CalculateCoveragePerBinsize.py
python CalculateCoveragePerBinsize.py hg38.100K.bin 100000 RUNX1.bedpe RUNX1.bedpe.bedgraph <br/>
对bed文件计算给定binsize的bedgraph结果<br/>

### SplitNLineBam.py
python SplitNLineBam.py input.bam 100000 outprefix <br/>
将给定的sam/bam文件拆分成多个bam文件，如果中间涉及了多位点比对，则文件必须按照name排序，然后有些文件会稍微多余给定的值<br/>

### downSample.py
python downSample.py -i ../7_reads_annotation/A549.CTCF.DNARNA.bedpe -n 15013571 -o A549.CTCF.DNARNA.bedpe -s 123456789 <br/>
对按行的文件取子样本<br/>

### fingBigData.sh
bash fingBigData.sh > user.text.size.txt <br/>
寻找大文件，代码中是寻找大于4G以上的文件，可以自定义 <br/>

### split_hisat2_bam.py
python split_hisat2_bam.py input.bam input <br/>
针对hisat2的比对结果拆分bam文件<br/>

### kmer.py
python kmer.py 100 10 full.fa genome.fastq.gz <br/>
将基因组文件以10bp为步长，提取100bp长度的序列片段，生成fastq文件<br/>

### codelibrary.java
java的常用工具集<br/>

### pair2single.py
python pair2single.py input.bam R1.bam R2.bam <br/>
将双端比对的结果拆分成单端比对的结果<br/>

### getRegionCoverage.py
python getRegionCoverage.py peak.file alignfile.bw peak.coverage <br/>
计算指定区域的平均测序深度<br/>