######### https://cloud.tencent.com/developer/article/2056780 # deseq2 差异基因分析，去除批次效应。
######### https://biohpc.cornell.edu/doc/RNA-Seq-2019-exercise3.html
######### https://blog.csdn.net/qazplm12_3/article/details/81221329
######### https://zhuanlan.zhihu.com/p/601284154
######### https://blog.csdn.net/zfyyzhys/article/details/140846587 # 画图好看
######### http://www.bioinformatics.com.cn/remove_batch_effects_by_sva_combat_t006 # 在线工具，画图好看
######### http://www.badd-cao.net/rank-in/submission.html

# 首先通过DESeq2和PCA分析识别批次效应，然后使用sva包的ComBat函数及DESeq2自身的调整来消除批次效应，最后探讨了limma包的removeBatchEffect方法。通过这些方法，可以确保后续差异表达分析的准确性，减少假阳性结果。

# 基因表达标准化，不同样品的测序量会有差异，最简单的标准化方式是计算
# counts per million (CPM)，即原始reads count除以总reads数乘以1,000,000。
# 这种计算方式的缺点是容易受到极高表达且在不同样品中存在差异表达的基因的影响；这些基因的打开或关闭会影响到细胞中总的分子数目，可能导致这些基因标准化之后就不存在表达差异了，而原本没有差异的基因标准化之后却有差异了。RPKM、FPKM和TPM是CPM按照基因或转录本长度归一化后的表达，也会受到这一影响。
# 为了解决这一问题，研究者提出了其它的标准化方法。量化因子 (size factor,SF)是由DESeq提出的。其方法是首先计算每个基因在所有样品中表达的几何平均值。每个细胞的量化因子(size factor)是所有基因与其在所有样品中的表达值的几何平均值的比值的中位数。由于几何平均值的使用，只有在所有样品中表达都不为0的基因才能用来计算。这一方法又被称为RLE (relative log expression)。
# 上四分位数 (upperquartile,UQ)是样品中所有基因的表达除以处于上四分位数的基因的表达值。同时为了保证表达水平的相对稳定，计算得到的上四分位数值要除以所有样品中上四分位数值的中位数（TCGA也提供了这个方法）。
# TMM是M-值的加权截尾均值。选定一个样品为参照，其它样品中基因的表达相对于参照样品中对应基因表达倍数的log2值定义为M-值。随后去除M-值中最高和最低的30%，剩下的M值计算加权平均值。每一个非参照样品的基因表达值都乘以计算出的TMM。这个方法假设大部分基因的表达是没有差异的。


calc_cpm <- function (expr_mat, spikes = NULL){
  norm_factor <- colSums(expr_mat[-spikes,])
  return(t(t(expr_mat)/norm_factor)) * 10^6
}

calc_sf <- function (expr_mat, spikes=NULL){
  geomeans <- exp(rowMeans(log(expr_mat[-spikes,])))
  SF <- function(cnts){
    median((cnts/geomeans)[(is.finite(geomeans) & geomeans >0)])
  }
  norm_factor <- apply(expr_mat[-spikes,],2,SF)
  return(t(t(expr_mat)/norm_factor))
}

calc_uq <- function (expr_mat, spikes=NULL){
  UQ <- function(x) {
    quantile(x[x>0],0.75)
  }
  uq <- unlist(apply(expr_mat[-spikes,],2,UQ))
  norm_factor <- uq/median(uq)
  return(t(t(expr_mat)/norm_factor))
}


# 校正批次效应的方法有很多，一项研究比较了6种去除批次效应的方法，其中包括ComBat方法（parametric prior method，ComBat_p和non-parametric method，ComBat_n）、代理变量法（Surrogate variable analysis，SVA）、基于比值的方法（Geometric ratio-based method，Ratio_G）、平均中心方法（Mean-centering，PAMR）和距离加权判别（Distance-weighted discrimination，DWD）方法，综合多种指标认为ComBat在精确性、准确性和整体性能方面（precision, accuracy and overall performance）总体优于其他方法。
# SVA的ComBat处理批次效应，SVA包中有两个函数可以用来校正批次效应：ComBat和ComBat_seq。ComBat使用参数或非参数经验贝叶斯模型，输入数据为干净的、标准化的表达数据，通常是芯片数据。ComBat_seq是一个使用负二项回归的ComBat改进模型，专门针对RNA-Seq count数据。我需要分析的数据便是RNAseq数据，后面例子便会详细介绍这个函数。
# raw counts矩阵（combat-seq函数处理）
# 芯片/FPKM等标准化矩阵（combat函数处理）

######### 标准化
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
sizefc = counts(dds, normalized=TRUE)


options(stringsAsFactors = F)
# BiocManager::install('airway')
# 加载airway数据集并转换为表达矩阵
library(airway,quietly = T)
library(gplots)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

data(airway)
class(airway)
rawcount <- assay(airway)
colnames(rawcount)
rawcount[1:4,1:4]
dim(rawcount)
group_list <- colData(airway)$dex
group_list
keep <- rowSums(rawcount>0) >= floor(0.75*ncol(rawcount))
table(keep)
filter_count <- rawcount[keep,]
filter_count[1:4,1:4]
dim(filter_count)
# 保存表达矩阵和分组结果
save(filter_count, group_list, file = "Step01-airwayData.Rdata") 


# 查看分组信息和表达矩阵数据
exprSet <- filter_count
dim(exprSet)
exprSet[1:6,1:6]
table(group_list)
exprSet[1:6,1:6]
table(group_list)
pheatmap(cor(exprSet))
batch=paste0('b',rep(1:2,each=4))
batch
dat=log2(edgeR::cpm(exprSet)+1)
ht_for_RNAcounts <- function(dat){
  dat[1:4,1:4] 
  colnames(dat)
  ac=data.frame(group=group_list,
                batch=batch)
  rownames(ac)=colnames(dat)
  pheatmap(cor(dat),annotation_col  =ac)
}
ht_for_RNAcounts(dat)

ct_with_batch = filter_count
x= sample(1:nrow(ct_with_batch))
y = which(batch=='b1')
ct_with_batch[x,y] = ct_with_batch[x,y]+100
dat=log2(edgeR::cpm(ct_with_batch)+1)
ht_for_RNAcounts(dat)


# 去除批次效应的差异分析
exprSet=ct_with_batch
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list, batch= batch )
dds <- DESeqDataSetFromMatrix(countData = exprSet, colData = colData, design = ~ group_list+batch )
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name=  "group_list_untrt_vs_trt")
summary(res)
save(res,file = 'group_list_untrt_vs_trt_rm_batch_deg.Rdata')

resOrdered <- res[order(res$padj),]
DEG =as.data.frame(resOrdered)
rm_batch_deg = na.omit(DEG)  
head(rm_batch_deg)
# 筛选上下调，设定阈值
fc_cutoff <- 1
fdr <- 0.05 
rm_batch_deg$regulated <- "normal" 
loc_up <- intersect(which(rm_batch_deg$log2FoldChange>log2(fc_cutoff)),
                    which(rm_batch_deg$padj<fdr))
loc_down <- intersect(which(rm_batch_deg$log2FoldChange< (-log2(fc_cutoff))),
                      which(rm_batch_deg$padj<fdr))
rm_batch_deg$regulated[loc_up] <- "up"
rm_batch_deg$regulated[loc_down] <- "down" 
table(rm_batch_deg$regulated)
rm_batch_deg$ENSEMBL = rownames(rm_batch_deg)
save(rm_batch_deg, file = "group_list_untrt_vs_trt_rm_batch_deg.result.Rdata")

# 保留批次效应的差异分析
exprSet=ct_with_batch
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet, colData = colData, design = ~ group_list)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name=  "group_list_untrt_vs_trt")
summary(res)
save(res,file = 'group_list_untrt_vs_trt_with_batch_deg.Rdata')

resOrdered <- res[order(res$padj),]
DEG =as.data.frame(resOrdered)
rm_batch_deg = na.omit(DEG)  
head(rm_batch_deg)
# 筛选上下调，设定阈值
fc_cutoff <- 1
fdr <- 0.05 
rm_batch_deg$regulated <- "normal" 
loc_up <- intersect(which(rm_batch_deg$log2FoldChange>log2(fc_cutoff)),
                    which(rm_batch_deg$padj<fdr))
loc_down <- intersect(which(rm_batch_deg$log2FoldChange< (-log2(fc_cutoff))),
                      which(rm_batch_deg$padj<fdr))
rm_batch_deg$regulated[loc_up] <- "up"
rm_batch_deg$regulated[loc_down] <- "down" 
table(rm_batch_deg$regulated)
rm_batch_deg$ENSEMBL = rownames(rm_batch_deg)
save(rm_batch_deg, file = "group_list_untrt_vs_trt_with_batch_deg.result.Rdata")


### 加载两个差异分析的结果
rm(list = ls())
load(file = "group_list_untrt_vs_trt_rm_batch_deg.result.Rdata")
rm_batch_deg1 = rm_batch_deg
load(file = "group_list_untrt_vs_trt_with_batch_deg.result.Rdata")
with_batch_deg1 = rm_batch_deg
tmp= merge(rm_batch_deg1,with_batch_deg1,by='ENSEMBL')
table(tmp$regulated.x,tmp$regulated.y)
balloonplot(table(tmp$regulated.x,tmp$regulated.y)) 




########## 展示如何去批次和去批次的效果，用limma包的removeBatchEffect函数来实现。
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("batch", "group_list"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=batch, shape = group_list)) +
 geom_point(size=3) +
 # xlim(-12, 12) +
 # ylim(-10, 10) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) +
 geom_text(aes(label=name),vjust=2)
# ggsave("myPCAWithBatchEffect.png")

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
pcaData <- plotPCA(vsd, intgroup=c("batch", "group_list"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=batch, shape = group_list)) +
 geom_point(size=3) +
 # xlim(-12, 12) +
 # ylim(-10, 10) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) +
 geom_text(aes(label=name),vjust=2)
# ggsave("myPCABatchEffectRemoved.png")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$batch, vsd$group_list, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
 clustering_distance_rows=sampleDists,
 clustering_distance_cols=sampleDists,
 col=colors)





















#####################################################
########   另一个案例                          #######
#####################################################
## 1. load数据
library(pasilla)
 
dataFilesDir = system.file("extdata", package = "pasilla", mustWork=TRUE)
pasillaSampleAnno = read.csv(file.path(dataFilesDir, "pasilla_sample_annotation.csv"))
count_df <- read.table(file.path(dataFilesDir,"pasilla_gene_counts.tsv"), sep="\t",header = TRUE, row.names=1)
count_matrix <- as.matrix(count_df)
 
#样品信息数据框，行名为表达矩阵的列名
colData = data.frame(
  condition = as.factor(c(rep("untreated",4),rep("treated",3))),
  type = as.factor(c("SR","SR","PE","PE","SR","PE","PE")))
rownames(colData) <- colnames(count_matrix)
 
## 2. 可视化查看批次效应
require(DESeq2)
library(ggplot2)
dds = DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design= ~ type + condition)
pcaData <- plotPCA(DESeqTransform(dds), intgroup=c("condition", "type"), returnData=TRUE)

# plotPCA {DESeq2}
percentVar <- round(100 * attr(pcaData, "percentVar"))
p1 = ggplot(pcaData, aes(PC1, PC2, color=condition, shape = type)) +
  geom_point(size=3) +
  # xlim(-4e+05, 5e+05) +
  # ylim(-4e+05, 2e+05) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)

## 3.使用sva包的ComBat 函数去除批次效应
library(sva)
count_mat <- ComBat(dat=count_matrix, batch=colData$type)
# ComBat:Adjust for batch effects using an empirical Bayes framework
count_mat[1:3,1:5]
table(count_mat < 0)
count_mat[count_mat<0] = 0

library(FactoMineR) # PCA函数
library(factoextra) # fviz_pca_ind函数
pre.pca <- PCA(t(count_mat),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= c("point", "text"),
             col.ind = colData$type,
             addEllipses = TRUE,
             legend.title="Group"  )




## 4. DESeq2包消除批次效应
# design 设计矩阵中加入引起批次效应的因素(SR or PE)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ type + condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="condition_untreated_vs_treated")
 
summary(res)
head(res)
resOdered <- res[order(res$padj),]
deg <- as.data.frame(resOdered)
#deg <- na.omit(deg)
dim(deg)
#write.csv(deg,file= "diff_deseq2.csv")
 
# 看看批次效应导致的差异表达基因,
# 这里是测序方法SR(single-end)和PE (paired-end)
batch_res <-  results(dds, name="type_SR_vs_PE")
summary(batch_res)
head(batch_res)
batchResOdered <- res[order(batch_res$padj),]
batchDeg <- as.data.frame(batchResOdered)
batchDeg <- na.omit(batchDeg)
dim(batchDeg)
head(batchDeg)


## 5. 用limma的removeBatchEffect去除批次效应
require(limma)
require(edgeR)  # DGEList来自edgeR包
expr_mat<- removeBatchEffect(count_matrix, batch=colData$type) 
table(expr_mat<0) #会出现负数！
expr_mat[which(expr_mat<0)]=0 # 不知是否合理！
expr_mat <- round(expr_mat) # read 为整数

# removeBatchEffect:This function is not intended to be used 
# prior to linear modelling. For linear modelling, 
# it is better to include the batch factors in the linear model.

# 查看去掉批次效应之后的数据
count_matrix[1:3,1:5]
expr_mat[1:3,1:5]
 
# 接着，查看数据PCA的分布是否消除了批次效应，
# 再用去除批次效应的expr_mat做下游分析
# dds = DESeqDataSetFromMatrix(countData = expr_mat,
#                              colData = colData,
#                              design= ~ condition +type)
# 
# pcaData <- plotPCA(DESeqTransform(dds), intgroup=c("condition", "type"), 
#                    returnData=TRUE)
# 
# # plotPCA {DESeq2}
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color=condition, shape = type)) +
#   geom_point(size=3) +
#   # xlim(-4e+05, 5e+05) +
#   # ylim(-4e+05, 2e+05) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   geom_text(aes(label=name),vjust=2)