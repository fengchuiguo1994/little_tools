# 可视化（浏览器）
### IGV是一个复合的支持所有常见类型的基因组数据的浏览器
##### IGV可视化二级结构配对关系
准备U1的基因组文件（前后延伸了5bp）<br/>
```
>U1
CAACATACTTACCTGGACGGGGTCGACGGCCGATCAAGAAGAGCCGTGGCCTAGGCCAATGGCCCACATTGCACTTGGTGGGCGCGTTGGCTTACCATCTCCCCAAGCGGGAGAGTGGACGTCATAATTTGTGGTAGAGGGGGTACGCGTTCGCGCGGCCCCTGCCAATT

实际U1序列
>U1
TACTTACCTGGACGGGGTCGACGGCCGATCAAGAAGAGCCGTGGCCTAGGCCAATGGCCCACATTGCACTTGGTGGGCGCGTTGGCTTACCATCTCCCCAAGCGGGAGAGTGGACGTCATAATTTGTGGTAGAGGGGGTACGCGTTCGCGCGGCCCCTGC
```
将实际的U1序列投入[RNAalifold webserver](http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAalifold.cgi)得到预测的RNA二级结构<br/>
```
UACUUACCUGGACGGGGUCGACGGCCGAUCAAGAAGAGCCGUGGCCUAGGCCAAUGGCCCACAUUGCACUUGGUGGGCGCGUUGGCUUACCAUCUCCCCAAGCGGGAGAGUGGACGUCAUAAUUUGUGGUAGAGGGGGUACGCGUUCGCGCGGCCCCUGC
...(((((..(((((((((.(((((...((.....))))))))))))(((((((((((((((..........)))))).))))))))).(((((((((.....)))))).))).)))).........))))).((((((.((((....))))))))))..  
```
得到的二级结构可以导入[R-chie](https://www.e-rna.org/r-chie/plot.cgi)生成模拟图。也可以用IGV实现<br/>
将二级结构转成点对点配对关系（python脚本）<br/>
```
aa = "U1"
bb = "...(((((..(((((((((.(((((...((.....))))))))))))(((((((((((((((..........)))))).))))))))).(((((((((.....)))))).))).)))).........))))).((((((.((((....)))))))))).."
# print(len(aa))
# print(len(bb))

mylist = []
for index,ii in enumerate(bb):
    if ii == "(":
        mylist.append(index+1)
    if ii == ")":
        p = mylist.pop()
        print("{0}\t{1}\t{2}".format(aa,p,index+1))
# print(len(mylist))

然后转化成实际位点：
awk '{if(NR==1){print $0}else{print $1"\t"5+$2-2"\t"5+$3-2}}' U1.bed > U1.cg.bed
```
![](figformanual\IGV.curve.png)<br/>
```
track graphType=arc
U1	33	39
U1	32	40
U1	28	41
U1	27	42
U1	26	43
U1	25	44
U1	24	45
U1	22	46
U1	21	47
U1	20	48
U1	19	49
U1	18	50
U1	65	76
U1	64	77
U1	63	78
U1	62	79
U1	61	80
U1	60	81
U1	59	83
U1	58	84
U1	57	85
U1	56	86
U1	55	87
U1	54	88
U1	53	89
U1	52	90
U1	51	91
U1	101	107
U1	100	108
U1	99	109
U1	98	110
U1	97	111
U1	96	112
U1	95	114
U1	94	115
U1	93	116
U1	17	118
U1	16	119
U1	15	120
U1	14	121
U1	11	131
U1	10	132
U1	9	133
U1	8	134
U1	7	135
U1	147	152
U1	146	153
U1	145	154
U1	144	155
U1	142	156
U1	141	157
U1	140	158
U1	139	159
U1	138	160
U1	137	161
```
##### IGV可视化bed文件
IGV支持bed3，bed4，bed6，bed12，narrowPeak（bed10）。其中介绍一下bed12<br/>
```
U1	28	143	27_Plus_A00988:31:HYVN7DSXY:3:1101:3703:16579	.	+	28	143	255,0,0	2	39,72	0,43
U1	14	138	96_Plus_A00988:31:HYVN7DSXY:3:1101:8666:20572	.	+	14	138	255,0,0	2	21,100	0,24
U1	27	138	161_Plus_A00988:31:HYVN7DSXY:3:1101:17065:7748	.	+	27	138	255,0,0	2	41,67	0,44
U1	73	165	205_Plus_A00988:31:HYVN7DSXY:3:1101:22218:16485	.	+	73	165	255,0,0	2	69,78	0,14
U1	17	144	460_Plus_A00988:31:HYVN7DSXY:3:1102:15998:34115	.	+	17	144	255,0,0	2	50,73	0,54
U1	20	134	563_Plus_A00988:31:HYVN7DSXY:3:1102:26549:27398	.	+	20	134	255,0,0	2	47,63	0,51
U1	33	105	567_Plus_A00988:31:HYVN7DSXY:3:1102:26928:28119	.	+	33	105	255,0,0	2	90,45	0,27
U1	27	139	593_Plus_A00988:31:HYVN7DSXY:3:1102:30463:16110	.	+	27	139	255,0,0	2	40,68	0,44
U1	26	133	687_Plus_A00988:31:HYVN7DSXY:3:1103:8223:29543	.	+	26	133	255,0,0	2	42,62	0,45
U1	17	138	719_Plus_A00988:31:HYVN7DSXY:3:1103:10728:28119	.	+	17	138	255,0,0	2	51,67	0,54
U1	20	138	791_Plus_A00988:31:HYVN7DSXY:3:1103:17047:16955	.	+	20	138	255,0,0	2	48,67	0,51
U1	27	145	843_Plus_A00988:31:HYVN7DSXY:3:1103:22761:10723	.	+	27	145	255,0,0	2	35,82	0,36
U1	22	138	955_Plus_A00988:31:HYVN7DSXY:3:1103:31530:36777	.	+	22	138	255,0,0	2	46,67	0,49
U1	27	144	987_Plus_A00988:31:HYVN7DSXY:3:1104:2275:19711	.	+	27	144	255,0,0	2	41,73	0,44
U1	21	146	1173_Plus_A00988:31:HYVN7DSXY:3:1104:19651:30608	.	+	21	146	255,0,0	2	47,75	0,50
U1	20	138	1286_Plus_A00988:31:HYVN7DSXY:3:1104:31801:21496	.	+	20	138	255,0,0	2	48,67	0,51
U1	23	143	1317_Plus_A00988:31:HYVN7DSXY:3:1105:3558:30138	.	+	23	143	255,0,0	2	43,72	0,48
U1	20	141	1456_Plus_A00988:31:HYVN7DSXY:3:1105:18548:6308	.	+	20	141	255,0,0	2	48,70	0,51
U1	17	146	1540_Plus_A00988:31:HYVN7DSXY:3:1105:26323:14199	.	+	17	146	255,0,0	2	50,75	0,54

chrom   start   end rid score   strand  start   end rgb_id  block_number    block_length    block_start
```
其中block_start，第一个一定是0，第二个是第二段的起点-第一段的起点。以此类推，都是n个block的起点减去第一个block的起点。

# 比对类工具
### 关于比对类工具处理兼并碱基
可以下载[水稻的细胞器基因组](https://rapdb.dna.affrc.go.jp/download/archive/Mt_Pt_genome.fasta)文件，里面包含了兼并碱基，也可以自定义。下载后，提取了一段reads做如下操作
```
@test1
ATTGGCGGGAGTATATTATGGCAGGATCAGTCACCTGGGCAAACCARCCC
+
ATTGGCGGGAGTATATTATGGCAGGATCAGTCACCTGGGCAAACCARCCC
@test2
ATTGGCGGGAGTATATTATGGCAGGATCAGTCACCTGGGCAAACCAACCC
+
ATTGGCGGGAGTATATTATGGCAGGATCAGTCACCTGGGCAAACCAACCC
@test3
ATTGGCGGGAGTATATTATGGCAGGATCAGTCACCTGGGCAAACCAGCCC
+
ATTGGCGGGAGTATATTATGGCAGGATCAGTCACCTGGGCAAACCAGCCC

bwa：
bwa index Mt_Pt_genome.fasta
bwa mem Mt_Pt_genome.fasta test.fq

hisat2:
hisat2-build Mt_Pt_genome.fasta Mt_Pt_genome
hisat2 -x Mt_Pt_genome -U test.fq

STAR:
mkdir -p starindex && STAR --runMode genomeGenerate --genomeDir starindex --genomeFastaFiles Mt_Pt_genome.fasta
STAR --genomeDir starindex/ --outFileNamePrefix startest --readFilesIn test.fq 
```
bwa的结果如图
![](figformanual\bwa.png)
将fastq文件中的兼并碱基变成了N。同时，基因组索引文件中，这个位置被记录成了G
hisat2的结果如图
![](figformanual\hisat2.png)
将fastq文件中的兼并碱基变成了A。而在基因组索引文件中，该位置被记录成N了。
STAR的结果如图
![](figformanual\STAR.png)
三个均无错配，说明可以识别兼并碱基
### blast
HAL1序列来源于RepBaseRepeatMaskerEdition-20181026.tar.gz，与[hg38 repeat](http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz)作比较。结果发现blast的比对结果短于注释的结果（bwa比对上的结果远远短于）。
```
makeblastdb -in hg38.fa -input_type fasta -dbtype nucl -parse_seqids -out hg38.blast.db -logfile hg38.log.txt
blastn -query HAL1.fa -out HAL1.blast.output -db hg38.blast.db -outfmt 6 -evalue 1e-5 -max_target_seqs 1000
```
输出文件说明
```
HAL1	chr12	74.464	979	172	43	1196	2153	49009079	49008158	1.40e-96	361
qacc：待比对的序列ID
sacc：比对上的ID
pident：相同匹配的百分比
length：比对上的长度
mismatch：错配数量
gapopen：缺口数量
qstart：待比对的序列起点
qend：待比对的序列终点
sstart：比对上的序列起点
send：比对上的序列终点
evalue：Expect value
bitscore：得分

# 亦可从下面选择任意自定义
qseqid qlen sseqid sgi slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname

blastn -query OTTP4.unmap.fa -out OTTP4.unmap.output -db hg38 -outfmt "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pidehg38 nidehg38 mismatch positive gapopen gaps ppos frames qframe sframe btop staxid scomname sblastname sskingdom staxids scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus ssciname sscinames" -evalue 1e-5 -max_target_seqs 2 -num_threads 12
```
### bwa
aln:
```
XT:U => Unique alignment
XT:R => Repeat
XT:N => Not mapped
XT:M => Mate-sw (成对比对的时候，当一端是被另一端救援而被确定比对位置的时候)
```
### STAR-Fusion 寻找融合RNA
```
STAR-Fusion --left_fq GM12878_R1.fq.gz --right_fq GM12878_R2.fq.gz --genome_lib_dir ~/Tools/STAR-Fusion/genome/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/ --CPU 20 --no_annotation_filter --examine_coding_effect --extract_fusion_reads --FusionInspector validate --denovo_reconstruct
```

# 定量工具
### featureCounts
[subread for featureCounts](https://github.com/DeskGen/subread).<br/>
可以对多位点比对的reads按百分数处理。单端数据：(针对RNA数据，因为一个reads只能有一个rna的来源，所以主要参数是--fraction -M)
```
featureCounts -a chr21.gtf -o chr21feature -t gene --largestOverlap --fraction -M -s 1  -J -G ~/Genome/GRCh3885/hg38.fa --donotsort -T 4 -R SAM --Rpath outtest/ --verbose -g gene_id chr21.bam

--fraction -M/-O ：平均分配多位点比对的结果
-g ：注意有些注释文件中可能不是gene_id，需要用-g修改

# 单端数据
featureCounts -a hg38.remove.gtf -o featureCounts.2.exon.txt -t exon -g gene_id -s 2 -R BAM -T 1 --Rpath exon2featureCounts --verbose uniqsplitfile/allrich-06.uniqmap.tmp.1.bam 

# 双端命令
featureCounts -a gencode_v26_Subread_to_DEXSeq.gtf -o SRR1797250.gene.count -t exon -g gene_id -p -T 1 SRR1797250.bam  
-p 将结果当成fragment来处理
-B 是不是双端都需要比对上(htseq默认是)
-P -d -D 设置插入片段大小
-C 不计算两端唯一不同染色体和位于同一染色体但链方向不一致的reads。
--donotsort 和htseq一样，都需要输入文件按name排序，不要相信提供的pos排序的参数。
--largestOverlap 分配给最长的比对的结果。

featureCounts -a ~/zhangyan/genome/hg38.analysisSet.chroms/hg38.remove.gtf -o featureCounts.${i}.${j}.txt -t ${j} -g gene_id -s ${i} -R BAM -T 1 --Rpath ${j}.${i}.featureCounts --verbose allrich-23.uniqmap.tmp.1.bam
featureCounts -a ~/zhangyan/genome/hg38.analysisSet.chroms/hg38.remove.gtf -o featureCounts.${i}.${j}.longest.txt -t ${j} -g gene_id --largestOverlap -s ${i} -R BAM -T 1 --Rpath ${j}.${i}.featureCounts.longest --verbose allrich-23.uniqmap.tmp.1.bam
```

# 序列拼接工具
### flash
[flash](http://ccb.jhu.edu/software/FLASH/index.shtml)
```
flash -M 150 -O --threads 12 --output-prefix out.flash --output-directory . read1.clean.fq read2.clean.fq
```
bug: 无法处理R1（R2）完全包含的reads。可能因为测序的偏差，导致R1和R2中的一段reads测序质量低被过滤掉。但是flash无法处理该情况。同时该情况在pandaseq中依然存在，而且pandaseq会输出一个明显错误的结果，flash输出为不拼接的结果。
```
R1: GTGAGCTGCCTTGGAAAAGGTTTGACATCATGGTCTCACCCTCCAGGCATTCGCAATGCTGTTGAAGCACTCTGGGCAATTCGGCTGGATTGCAACAGCCTCCTCGTTCTTCGCGATGCACATGTCAAACTCTCGTAGCTAAACCAAATC
R2: 
ATTTGGTTTAGCTACGAGAGGTTGACATGTGCATCGCGAAGAACGAGGAGGCTGTTGCAATCCAGCCGAATTGCCCAGAGTGCTTCAACAGCATTGCGAATGCCTGGAGGGTGAGACCATGATGTCAAACCTTTTCCAAGGC

R1: GTGAGCTGCCTTGGAAAAGGTTTGACATCATGGTCTCACCCTCCAGGCATTCGCAATGCTGTTGAAGCACTCTGGGCAATTCGGCTGGATTGCAACAGCCTCCTCGTTCTTCGCGATGCACATGTCAACCTCTCGTAGCTAAACCAAATC
R2: 
ATTTGGTTTAGCTACGAGAGGTTGACATGTGCATCGCGAAGAACGAGGAGGCTGTTGCAATCCAGCCGAATTGCCCAGAGTGCTTCAACAGCATTGCGAATGCCTGGAGGGTGAGACCATGATGTCAAACCTTTTCCAAGGC
```

# k-mer
### 得到kmer序列
[dsk](http://minia.genouest.org/dsk/)
```
i=sample
dsk  -file ${i}.1.fq,${i}.2.fq -kmer-size 20 -out $i.20.h5 -out-tmp $i
dsk2ascii -file $i.20.h5 -out $i.20.txt
```

# 变异
### 变异注释
[网页学习链接Annovar，SnpEff，Oncotator](http://yangli.name/2016/05/15/20160515annotation/)

# 注释
### 转录组注释
[InterProscan](http://www.ebi.ac.uk/interpro/download/)，[GitHub](https://github.com/ebi-pf-team/interproscan)
```
module load interproscan/5.48-83.0 && interproscan.sh --appl Pfam -t n -dp -i Pfam/dir${i}/part.${i}.fa --iprlookup -f tsv -o Pfam/dir${i}/part.${i}.fa.interproscan -cpu 5 -goterms -pa
可以移除参数--appl Pfam，对所有的数据库进行注释，默认是对蛋白序列，-t n修改为DNA/RNA核苷酸序列
```

# Linux
### sort （Linux常用命令）
看起来sort是忽略了字母大小写排序，目前还没仔细研究，日后再仔细研究。现在碰到sort请全部用LANG=C sort来代替。运行速度更快，跟其它软件sort之后的结果一致。
```
echo -e "snoR1\nSNOR75\nsnoR1\nsnoR75" > test
sort test 
    snoR1
    snoR1
    snoR75
    SNOR75

LANG=C sort test
    SNOR75
    snoR1
    snoR1
    snoR75
```

# 统计学
### 变异系数（CV）
1：标准差与平均数的比值称为变异系数，记为C·V。
2：变异系数可以消除单位和（或）平均数不同对两个或多个资料变异程度比较的影响。
```
aa = c(5.8,4.6,4.9)
sd(aa)/mean(aa)
```

# RNA二级结构
### [RNAup](http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAup.cgi)
评估两个RNA的结合能
### [RNAcofold](http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAcofold.cgi)
给出两个RNA的结合位点和方式
```
>U4
ATCTTTGCGCTTGGGGCAATGACGCAGCTTATGAGGTTATACCGAGGCGCGTCAATTGCTGGTTGAAAACTATTTCCAAACCCCCTCTTTGGCCCACCCGGGCCACTGAGAATTTCTGGAAGGGCTCTCTCTCTTCGGAGAGAGTAAAGC
>U6
TCTCTTCGGAGACATCCGATAAAATTGGAACGATACAGAGAAGATTAGCATGGCCCCTGCGCAAGGATGACACGCACAAATCGAGAAATGGTCCAAATTTTT
```
得到
```
>U46
AUCUUUGCGCUUGGGGCAAUGACGCAGCUUAUGAGGUUAUACCGAGGCGCGUCAAUUGCUGGUUGAAAACUAUUUCCAAACCCCCUCUUUGGCCCACCCGGGCCACUGAGAAUUUCUGGAAGGGCUCUCUCUCUUCGGAGAGAGUAAAGC&UCUCUUCGGAGACAUCCGAUAAAAUUGGAACGAUACAGAGAAGAUUAGCAUGGCCCCUGCGCAAGGAUGACACGCACAAAUCGAGAAAUGGUCCAAAUUUUU
((((((((((..(((((..((((((.(((.....(((...)))..)))))))))..((((((((........((((((......(((..((((((....))))))..)))......))))))...(((((.((((((((((((.......&))))))))))))..((((((...)))))).......))))).))))))))..))))).))))))))))...............((((((.......))))))  
```
导入[R-chie](https://e-rna.org/r-chie/plot.cgi)得到
![](figformanual\U46.png)
