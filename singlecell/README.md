## 通用
#### ChIAPIPE数据分析需要识别正确的测序标签，从geo上下载的数据已经丢失read id，变为SRR的序列，人为构建正确的read id。
```
python /data/home/ruanlab/huangxingyu/Tools/littletools/sra2readsid.py -r1 SRR17666401_1.fastq.gz -r2 SRR17666401_2.fastq.gz -o1 SRR17666401.name_1.fastq -o2 SRR17666401.name_2.fastq && pigz -p 10 SRR17666401.name_1.fastq && pigz -p 10 SRR17666401.name_2.fastq
```

## scHiC
#### [BandNorm](https://sshen82.github.io/BandNorm/index.html)<br/>
对scHiC数据聚类
```
python BandNormbedpe2mat.py file.bedpe file.mat

file.bedpe数据格式示例
chr1	41526218	41526368	chr1	42735421	42735468	SCG0092_AAACAGCCAGACAAAC-1	101302810:22096:8454144:58786689:1433:15646:21324	37	+	-
chr3	145107632	145107658	chr3	145118537	145118619	SCG0092_AAACAGCCAGACAAAC-1	103933621:24154:8454274:58852226:1440:31584:15859	37	-	+
chr6	29590106	29590125	chr6	29590171	29590226	SCG0092_AAACAGCCAGACAAAC-1	107906240:24154:8519809:58786689:1450:30219:27806	20	+	-
chr10	68252340	68252481	chr10	68550776	68550830	SCG0092_AAACAGCCAGACAAAC-1	10895412:22098:8454274:58852226:1130:7699:15389	37	-	+
chr8	84650218	84650372	chr8	84650413	84650437	SCG0092_AAACAGCCAGACAAAC-1	11440121:22096:8519680:58852226:1131:23393:30859	25	+	-
chr16	22429590	22429658	chr16	22615984	22616057	SCG0092_AAACAGCCAGACAAAC-1	115511568:24154:8519809:58786689:1470:14091:21574	37	+	-
chr7	68034221	68034244	chr7	68034294	68034421	SCG0092_AAACAGCCAGACAAAC-1	117345993:22098:8454274:58852226:1475:10700:14137	37	+	-
chr2	61743478	61743506	chr2	61749851	61749989	SCG0092_AAACAGCCAGACAAAC-1	119689623:22096:8454144:58786689:1503:18005:21872	37	+	+
chr8	3395221	3395259	chr8	3395726	3395764	SCG0092_AAACAGCCAGACAAAC-1	120804683:24154:8519809:58786689:1506:16631:21183	37	-	-
chr19	34254557	34254615	chr19	34254890	34255027	SCG0092_AAACAGCCAGACAAAC-1	121506217:22098:8454274:58852226:1508:26332:16282	37	+	+
```

## scRNA-Seq
```
library(Signac)
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(clustree)
library(patchwork)

#加载数据
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")

# 创建Seurat对象
seurat <- CreateSeuratObject(counts, project = "SCG0164", min.cells = 5)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)

#降维、聚类
seurat <- RunPCA(seurat)
seurat <- FindNeighbors(seurat, dims = 1:20)
# seurat <- FindClusters(seurat, resolution = 0.6)
seurat <- FindClusters(seurat, resolution = c(seq(0,2,.1)))
seurat <- RunUMAP(seurat,dims = 1:20)

clustree(seurat@meta.data, prefix = "RNA_snn_res.", node_colour_aggr = "median")
# clustree(seurat, prefix = "RNA_snn_res.", node_colour = "purple", node_size = 10, node_alpha = 0.8)

# 可视化结果
DimPlot(seurat, reduction = "umap", label = TRUE, group.by = "RNA_snn_res.0.4") + DimPlot(seurat, reduction = "umap", label = TRUE, group.by = "RNA_snn_res.0.8")

seurat<- AddMetaData(seurat,seurat@reductions$umap@cell.embeddings, col.name = c("UMAP_1","UMAP_2"))
seurat<- AddMetaData(seurat,seurat@reductions$pca@cell.embeddings, col.name = colnames(seurat@reductions$pca@cell.embeddings))
clustree_overlay(seurat, prefix = "RNA_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2")
clustree_overlay(seurat, prefix = "RNA_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2", use_colour = "points", alt_colour = "blue")
clustree_overlay(seurat, prefix = "RNA_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2", use_colour = "points", alt_colour = "blue", label_nodes = TRUE,plot_sides = TRUE)

clustree(seurat, prefix = "RNA_snn_res.") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "blue", high = "red")+
  theme(legend.position = "bottom")

clustree(seurat) + 
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set3")


sapply(grep( "^RNA_snn_res",colnames(seurat@meta.data),value = TRUE), function( x) length(unique(seurat@meta.data[, x])))
clip()


# 整合数据
ifnb.list = list(SCG0054=SCG0054flt,SCG0074=SCG0074flt,SCG0075=SCG0075flt,SCG0076=SCG0076flt,SCG0077=SCG0077flt,SCG0078=SCG0078flt,SCG0079=SCG0079flt,SCG0080=SCG0080flt)
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
```

## 时序分析
```
#####################################################################################
##0. 包的安装和加载
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("org.Rn.eg.db")
#BiocManager::install("org.At.tair.db")
#BiocManager::install("GO.db")
#BiocManager::install("monocle")
#BiocManager::install("AUCell")
#install.packages("rjson")
#install.packages("stringr")
#载入R包；
#devtools::load_all("D:/Program_Files/R-4.1.1/library/monocle")
library(monocle)
library(Seurat)
library(AUCell)
library(patchwork)
library(ggplot2)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.At.tair.db)
library(dplyr)
library(GO.db)
library(rjson)
library(stringr)

####################################################################################
##1. 数据加载、提取、构建CDS对象
##1.1 从seurat对象中提取部分细胞作为拟时分析的输入数据
#加载seurat对象的数据
setwd("C:/Genedenovo/单细胞培训班/单细胞培训班/单细胞三天培训班/NO.3 单细胞转录组高级分析及绘图/8.monocle实操脚本及数据")
#load("obj.Rda")
load("harmony_obj.Rda")
seurat_object@meta.data$celltype = as.vector(seurat_object@meta.data$seurat_clusters)
seurat_object@meta.data$celltype[ seurat_object@meta.data$seurat_clusters == "0"] = "Naive CD4 T"
seurat_object@meta.data$celltype[ seurat_object@meta.data$seurat_clusters == "1"] = "Memory CD4 T"
seurat_object@meta.data$celltype[ seurat_object@meta.data$seurat_clusters == "4"] = "CD8 T"
seurat_object@meta.data$celltype = factor(seurat_object@meta.data$celltype)
#统计细胞类型在样本中的分布
Idents(seurat_object) = "celltype"
table(seurat_object$orig.ident,seurat_object$celltype)
#根据细胞数量提取特定的细胞类型进行拟时分析
Tcell <- subset(seurat_object,idents = c("Naive CD4 T","Memory CD4 T","CD8 T"))
dim(Tcell)
Tcell=subset(Tcell,orig.ident=="stimulus")
dim(Tcell)
#提取亚群的umap图分布
DimPlot(seurat_object, reduction = "umap") + DimPlot(Tcell, reduction = "umap")

#为了减少内存占用，可以保存提取数据的seurat对象，并且删除原来的seurat对象
save(Tcell,file = "Tcell.Rda")
rm(seurat_object)

#从亚群seurat对象中提取拟时分析所需的数据
exp <- Tcell@assays[["RNA"]]@counts
fdata <- data.frame(gene_short_name = row.names(Tcell), row.names = row.names(Tcell))
pdata <- Tcell@meta.data

##1.2 CDS对象的创建
#将基因特征文件和细胞表型文件重新写一个对象
fd <- new("AnnotatedDataFrame", data = fdata)
pd <- new("AnnotatedDataFrame", data = pdata)

#构建monocle分析的CDS对象
CDS <- newCellDataSet(cellData = exp,phenoData = pd,featureData = fd)

#计算size factors 和 dispersions（离差），用于后期分析；
#结果：在phenoData表格添加1列Size_Factor；
CDS <- estimateSizeFactors(CDS)
CDS <- estimateDispersions(CDS)

#fData()函数用于提取CDS对象中的基因注释表格，得到的结果为数据框；
#pData()函数作用类似，提取CDS对象中的细胞表型表格；
head(pData(CDS))
head(fData(CDS))
head(dispersionTable(CDS))

#保存创建的CDS对象
save(CDS,file = "CDS_new.Rdata")
rm(list = ls())
load("CDS_new.Rdata")
dim(CDS)

############################################################################################
##2. 差异分析寻找高变基因
#detectGenes()函数：同时统计表达当前基因的细胞数量和细胞表达的基因数；
#min_expr参数用于设置检测阈值，比如min_expr = 0.1表示当前基因的表达量超过0.1才会纳入统计；
CDS <- detectGenes(CDS, min_expr = 0.1)
#过滤低表达的基因，以降低噪音和减少差异分析计算量
expressed_genes <-  row.names(subset(fData(CDS),num_cells_expressed >= 20));length(expressed_genes)
#不同细胞类型的差异分析
clustering_DEG_genes <-differentialGeneTest(CDS[expressed_genes,],fullModelFormulaStr = '~celltype')

############################################################################################
##3. 时间轨迹及其差异分析
##3.1 构建细胞轨迹
#第一步：选择用于构建细胞轨迹的基因集；
#选择Top2000差异基因作为排序基因；
ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:2000]
#fData(CDS)$use_for_ordering <-fData(CDS)$num_cells_expressed > 0.1 * ncol(CDS)
#将差异基因储存到CDS对象中；
CDS <- setOrderingFilter(CDS,ordering_genes = ordering_genes)
plot_ordering_genes(CDS)

#第二步: 数据降维
#降维函数与上文的t-SNE一致，但降维算法这里用的是“DDRTree”，
CDS <- reduceDimension(CDS, method = 'DDRTree')

#第三步: 构建细胞分化轨迹
#按照分化轨迹排序细胞；
CDS <-orderCells(CDS)

#绘制细胞分化轨迹：
#按“Pseudotime”分组；
plot_cell_trajectory(CDS, color_by = "Pseudotime")
#按“State”分组；
plot_cell_trajectory(CDS, color_by = "State")
#按seurat分群结果分组
plot_cell_trajectory(CDS, color_by = "seurat_clusters")
#按细胞类型分组
plot_cell_trajectory(CDS, color_by = "celltype")

# save(CDS,file = "CDS_pseudotime.Rda")
load("CDS_pseudotime.Rda")
##3.2 比较细胞分化轨迹进程中功能基因的表达差异
#主要用到sm.ns()函数根据表达量拟合曲线；
diff_test_res <- differentialGeneTest(CDS[expressed_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)")

head(diff_test_res[,c("gene_short_name", "pval", "qval")])
#按q值从小到大排序后，查看最显著的前4个基因的拟时间趋势
sig_gene_names1 <- rownames(diff_test_res[order(diff_test_res$qval)[1:4],])
plot_genes_in_pseudotime(CDS[sig_gene_names1,], color_by = 'Pseudotime')
plot_genes_in_pseudotime(CDS[sig_gene_names1,], color_by = 'State')
#也可以查看比较关注基因的拟时间趋势
sig_gene_names1 <- c("GZMB","CXCL10","RPL13","RPL7")
plot_genes_in_pseudotime(CDS[sig_gene_names1,], color_by = 'Pseudotime')

##3.3 对比较关注的基因使用ggplot2重新绘制非线性拟合曲线，可用于文章发表
#从CDS对象中提取数据
data <- data.frame(pData(CDS)[,c("Pseudotime","State")],t(as.matrix(CDS[c("GZMB","CXCL10","RPL13","RPL7"),])))
#将宽数据整理成长数据，方便使用ggplot2绘制多曲线图
data1=data.frame()
for (i in 1:(ncol(data)-2)){
  a=data.frame(cell = rownames(data),Pseudotime=data[,1],Expression=data[,2+i],gene=rep(colnames(data)[2+i],nrow(data)))
  data1=rbind(data1,a)
} 
#按照指定顺序绘制不同gene曲线
data1$gene <- factor(data1$gene,levels = unique(data1$gene),ordered = T)
#数据映射到图形上
p1 <- ggplot(data1)+geom_smooth(aes(Pseudotime,Expression,color=gene),method = "loess",se = F,size = 1.2);p1
p2 <- p1+scale_color_manual(values = c("red","OrangeRed","green","SeaGreen"));p2
p3 <- p2+labs(color="Gene Name",x="Pseudo-time",y="Gene Expression",
              title="Pseudo-time curve of different genes");p3
mytheme<-theme_bw()+theme(text = element_text(family = "sans"),  #调整图形字体型号
                          plot.title = element_text(size = rel(1.2),hjust = 0.5),  #图形标题字体大小及居中
                          axis.title = element_text(size = rel(1)),  #坐标轴标题字体大小
                          axis.text = element_text(size=rel(0.8),colour = "black"),  #坐标轴刻度的字体大小
                          legend.text = element_text(size = rel(0.8)),  #图例标题字体大小
                          legend.title = element_text(size = rel(1)),  #图例内容字体大小
                          plot.margin=unit(x=c(0.2,0.2,0.2,0.2),units="inches"),  #图形外边框的间距
                          panel.grid = element_blank()) #去除网格线
p4<-p3+mytheme;p4
ggsave(p4,filename = "curve_gene.pdf",width = 15,height = 15,units = "cm")

##3.4 差异基因的拟时表达模式聚类分析
#提取差异基因；
sig_gene_names2 <- row.names(subset(diff_test_res, qval < 0.0001))
#绘制拟时间差异基因表达谱热图；
plot_pseudotime_heatmap(CDS[sig_gene_names2[1:100],],num_clusters = 3,
                        show_rownames = F)

######################################################################################################
##4.单细胞轨迹分支分析 
#当细胞分化轨迹出现分支的时候，意味着细胞将面临不同的分化命运“fate”，接下来主要分析分支事件，比如沿着分化轨迹，基因的表达量如何变化？
#不同分支之间的差异基因有哪些？

#Monocle 提供一个特定的统计检验方法: branched expression analysis modeling（BEAM）.
##4.1 BEAM检验
#使用BEAM()函数对基因进行注释；
#BEAM函数的输入对象： 完成拟时间分析的CellDataSet且轨迹中包含1个分支点；
#返回一个包含每个基因significance scores 的表格,若得分是显著的则表明该基因的表达是与分支相关的（branch-dependent）。

BEAM_res <- BEAM(CDS[expressed_genes,], branch_point = 1, progenitor_method = "duplicate")
#按照qval升序排列；
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
head(BEAM_res)

##4.2轨迹分支表达分析
#使用pheatmap包绘制分支表达量热图；
plot_genes_branched_heatmap(CDS[row.names(subset(BEAM_res,qval < 1e-5)),],branch_point = 1,num_clusters = 4
                            ,use_gene_short_name = F,show_rownames = F)

#使用 plot_genes_branched_pseudotime() 函数绘制拟合曲线
plot_genes_branched_pseudotime(CDS[rownames(BEAM_res)[1:4],],branch_point = 1,color_by = "Pseudotime",ncol = 1)

###################################################################################################
##5. 基因集分析
##5.1 免疫相关的GO term基因集的查找
#获取人类的所有GO term基因集
source("getGoTerm.R")
GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
save(GO_DATA, file = "GO_DATA.Rda")

#根据关键词查找相关通路的ID
load("GO_DATA.Rda") # 载入数据 GO_DATA
findGO("T cell") # 寻找含有指定关键字的 pathway name 的 pathway
findGO("CXCL10", method = "gene") # 寻找含有指定基因名的 pathway

#T cell differentiation involved in immune response
T_cell_differentiation <- getGO("GO:0002292");T_cell_differentiation

#T cell receptor signaling pathway
T_cell_costimulation <- getGO("GO:0031295");T_cell_costimulation

#T cell lineage commitment
T_cell_lineage_commitment <- getGO("GO:0002360");T_cell_lineage_commitment


##5.2 免疫相关的KEGG通路基因集  https://www.kegg.jp/kegg/pathway.html  Search=T cell immune
#Systemic lupus erythematosus
download.file("http://togows.dbcls.jp/entry/pathway/hsa05322/genes.json", "hsa05322.json")
Systemic_lupus_erythematosus = fromJSON(file ="hsa05322.json")
SLE_geneset = list(as.character(sapply(Systemic_lupus_erythematosus[[1]], function(x) sapply(strsplit(x[1], ";"), function(x) x[1]))))
SLE_geneset

#T cell receptor signaling pathway
download.file("http://togows.dbcls.jp/entry/pathway/hsa04660/genes.json", "hsa04660.json")
T_cell_receptor_signaling_pathway = fromJSON(file ="hsa04660.json")
TCR_geneset = list(as.character(sapply(T_cell_receptor_signaling_pathway[[1]], function(x) sapply(strsplit(x[1], ";"), function(x) x[1]))))
TCR_geneset

#Th17 cell differentiation
download.file("http://togows.dbcls.jp/entry/pathway/hsa04659/genes.json", "hsa04659.json")
Th17_cell_differentiation = fromJSON(file ="hsa04659.json")
TCD_geneset = list(as.character(sapply(Th17_cell_differentiation[[1]], function(x) sapply(strsplit(x[1], ";"), function(x) x[1]))))
TCD_geneset

#方法二
#source("hsa_kegg_path_db.R")
#findKEGG("T cell")
#getKEGG("hsa04660")

#MsigDB通路搜索
source("hsa_MsigDB_path_db.R")
findMsigDB("signaling")
Il6_Jak_Stat3_Signaling <- getMsigDB("Il6 Jak Stat3 Signaling");Il6_Jak_Stat3_Signaling

##5.3 基因集score值计算
load("Tcell.Rda")
#计算每个细胞的基因集score值。具体过程为先计算目标基因集在每个细胞中所有目标基因的平均值，再根据平均值把表达矩阵切割成若干个文件，
#然后从切割后的每一个文件中随机抽取对照基因作为背景值（基因集外的基因，默认取100个）。
#最后所有的目标基因算一个平均值，所有的背景基因算一个平均值，两者相减就是改基因集在某一个细胞的score值。
#所以会出现相同基因集在相同的细胞打分，也可能会产生不同的score值。
#免疫相关的GO基因集评分
Tcell <- AddModuleScore(Tcell,features = T_cell_differentiation,name = "T_cell_differentiation")
Tcell <- AddModuleScore(Tcell,features = T_cell_costimulation,name = "T_cell_costimulation")
Tcell <- AddModuleScore(Tcell,features = T_cell_lineage_commitment,name = "T_cell_lineage_commitment")

#免疫相关的KEGG基因集评分
Tcell <- AddModuleScore(Tcell,features = SLE_geneset,name = "Systemic_lupus_erythematosus")
Tcell <- AddModuleScore(Tcell,features = TCR_geneset,name = "T_cell_receptor_signaling_pathway")
Tcell <- AddModuleScore(Tcell,features = TCD_geneset,name = "Th17_cell_differentiation")

Tcell <- AddModuleScore(Tcell,features = Il6_Jak_Stat3_Signaling,name = "Il6_Jak_Stat3_Signaling")

##5.4 不同基因集的拟时间曲线图
cell_gene_score <- FetchData(Tcell,vars = c(colnames(Tcell@meta.data)[-7:-1]))
cell_gene_score1 = cell_gene_score[rownames(pData(CDS)),]
data <- data.frame(pData(CDS)[,c("Pseudotime","State")],cell_gene_score1)
data <- data.frame(pData(CDS)[,c("Pseudotime","State")],cell_gene_score)
colnames(data)[-2:-1] <- sub("1$","",str_replace_all(colnames(data)[-2:-1],"_", " "))

#将宽数据整理成长数据，方便使用ggplot2绘制多曲线图
data1=data.frame()
for (i in 1:(ncol(data)-2)){
  a=data.frame(cell = rownames(data),Pseudotime=data[,1],Expression=data[,2+i],gene=rep(colnames(data)[2+i],nrow(data)))
  data1=rbind(data1,a)
} 
#按照指定顺序绘制不同gene set曲线
data1$gene <- factor(data1$gene,levels = unique(data1$gene),ordered = T)
#数据映射到图形上
p1 <- ggplot(data1)+geom_smooth(aes(Pseudotime,Expression,color=gene),method = "loess",se = F,size = 1.2);p1
p2 <- p1+scale_color_manual(values = c("red","OrangeRed","green","SeaGreen","#FFD700","#1E90FF","blue"));p2
p3 <- p2+labs(color="GO Term/\nKEGG Pathway",x="Pseudo-time",y="Gene set score",
              title="Pseudo-time curve of different gene set score");p3
mytheme<-theme_bw()+theme(text = element_text(family = "sans"),  #调整图形字体型号
                          plot.title = element_text(size = rel(1.2),hjust = 0.2),  #图形标题字体大小及居中
                          axis.title = element_text(size = rel(1)),  #坐标轴标题字体大小
                          axis.text = element_text(size=rel(0.8),colour = "black"),  #坐标轴刻度的字体大小
                          legend.text = element_text(size = rel(0.8)),  #图例标题字体大小
                          legend.title = element_text(size = rel(1)),  #图例内容字体大小
                          plot.margin=unit(x=c(0.2,0.2,0.2,0.2),units="inches"),  #图形外边框的间距
                          panel.grid = element_blank()) #去除网格线
p4<-p3+mytheme;p4
ggsave(p4,filename = "curve_gene_set_score.pdf",width = 18,height = 15,units = "cm")

##5.5 小提琴图展示
VlnPlot(Tcell,features = "T_cell_differentiation1",
        pt.size = 0.1, adjust = 1,group.by = "celltype")

#5.6 基因得分的umap映射图
mydata <- FetchData(Tcell,vars = c("UMAP_1","UMAP_2","T_cell_differentiation1"))
p1 <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = T_cell_differentiation1))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
  colours = c("#333366","#6666FF","#CC3333","#FFCC33"));p1
p2 <- p1 + labs(colour="Gene set score",title = "UMAP graph of T cell differentiation");p2
p3 <- p2 + theme_classic() + theme(text = element_text(family = "sans"),  #调整图形字体型号
                              plot.title = element_text(size = rel(1.4),hjust = 0.3,face = "bold"),  #图形标题字体大小及居中
                              axis.title = element_text(size = rel(1.3),colour = "black"),  #坐标轴标题字体大小
                              axis.text = element_text(size=rel(1.1),colour = "black"),  #坐标轴刻度的字体大小
                              legend.text = element_text(size = rel(1)),  #图例标题字体大小
                              legend.title = element_text(size = rel(1.1)),  #图例内容字体大小
                              plot.margin=unit(x=c(0.2,0.2,0.2,0.2),units="inches"),  #图形外边框的间距
                              panel.grid = element_blank());p3 #去除网格线

p4 <- p3+DimPlot(Tcell,reduction = "umap",group.by = "celltype",label = T,pt.size = 1.2,label.size = 3);p4
ggsave(p4,filename = "UMAP_gene_score.pdf",width = 25,height = 12,units = "cm")

##5.7 箱型散点图
data<- FetchData(Tcell,vars = c("celltype","T_cell_receptor_signaling_pathway1"))
p1 <- ggplot(data, aes(x=celltype,y=T_cell_receptor_signaling_pathway1)) + geom_jitter(col="#00000033", pch=19,cex=2, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),aes(color = factor(celltype)));p1
p2 <- p1 + labs(x=NULL,y="Gene set score",title = "T cell receptor signaling pathway",color = "Cell Type");p2
p3 <-  p2 + theme_bw() + theme(text = element_text(family = "sans"),  #调整图形字体型号
                               plot.title = element_text(size = rel(1.2),hjust = 0.5),  #图形标题字体大小及居中
                               axis.title = element_text(size = rel(1)),  #坐标轴标题字体大小
                               axis.text = element_text(size=rel(0.8),colour = "black"),  #坐标轴刻度的字体大小
                               axis.text.x = element_text(size = rel(1.1)),
                               legend.text = element_text(size = rel(0.8)),  #图例标题字体大小
                               legend.title = element_text(size = rel(1)),  #图例内容字体大小
                               plot.margin=unit(x=c(0.2,0.2,0.2,0.2),units="inches"),  #图形外边框的间距
                               panel.grid = element_blank());p3 #去除网格线
ggsave(p3,filename = "boxplot_gene_score.pdf",width = 15,height = 15,units = "cm")

######################################################################################
##6 SCENIC计算每个细胞的TF regulons的AUC值（感兴趣的老师可以尝试）（不是本节课的重点）
##该部分消耗内存和CPU极大，有32G内存电脑的老师可以尝试运行，运行时间大约2小时左右。输出结果就是文件中的int和output两个文件，
##6.1 包的安装与加载
#BiocManager::install(version = "3.15")
#BiocManager::install(c("AUCell", "RcisTarget"))
#BiocManager::install(c("GENIE3"))
#BiocManager::install(c("zoo", "mixtools", "rbokeh"))
#BiocManager::install(c("DT", "NMF", "pheatmap", "R2HTML", "Rtsne", "doRNG"))
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
#devtools::install_github("aertslab/SCENIC")
#install.packages("dplyr")
#install.packages("doMC", repos="http://R-Forge.R-project.org")

#library(SCENIC)
#library(dplyr)
#library(Seurat)
#library(foreach)

##6.2 下载人类的SCENIC数据库
#下载500bp upstream
#RcisTarget_hgnc_500bp <- "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather"
#download.file(RcisTarget_hgnc_500bp, destfile=basename(RcisTarget_hgnc_500bp))
#下载TSS+/-10kbp
#RcisTarget_hgnc_10kbp <- "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather"
#download.file(RcisTarget_hgnc_10kbp, destfile=basename(RcisTarget_hgnc_10kbp))

##6.3 从seurat对象中提取数据用于SCENIC的输入
#load("Tcell.Rda")
#exprMat <- as.matrix(Tcell@assays$RNA@counts)
#dim(exprMat)
#cell.meta <- data.frame(Tcell@meta.data)
#colnames(cell.meta)[which(colnames(cell.meta)=="orig.ident")] <- "sample"
#colnames(cell.meta)[which(colnames(cell.meta)=="seurat_clusters")] <- "cluster"
#cell.meta <- cell.meta[,c("sample","cluster")]

##6.4 初始设置
#scenicOptions <- initializeScenic(org = "hgnc",dbDir = ".",nCores = 4)
#saveRDS(cell.meta,file = "int/cell.meta.Rds")
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cell.meta.Rds"

##6.5 构建网络
#genesKept <- geneFiltering(exprMat = exprMat,scenicOptions = scenicOptions,minCountsPerGene = 3 * 0.01 *ncol(exprMat),minSamples = ncol(exprMat)*0.01)
#exprMat_filter <- exprMat[genesKept,]
#runCorrelation(exprMat_filtered = exprMat_filter,scenicOptions = scenicOptions)
#exprMat_filtered_log <- log2(exprMat_filter+1)
#runGenie3(exprMat = exprMat_filtered_log,scenicOptions = scenicOptions)

##6.6 基因调控网络的构建和评分
#exprMat_log <- log2(exprMat+1)
#scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

########################################################################
rm(list = ls())
##7. TF regulons的AUC值拟时间曲线
#基于样本中的regulons基因表达排名（gene expression rank），使用曲线下面积来评估输入基因集是否在样本的前5%表达基因内富集
##7.1 整理数据
#加载SCENIC计算的regulonAUC的rds文件
regulonAUC=readRDS("3.4_regulonAUC.Rds")
#从rsd文件种提取每个细胞的不同regulons的AUC值
AUC_cell_data=t(regulonAUC@assays@data@listData[["AUC"]])
#加载CDS对象
load("CDS_pseudotime.Rda")
head(pData(CDS)) #查看CDS对象的细胞表型
#将每个细胞的AUC值与拟时间值整理成一个数据框
AUC_cell_data1=data.frame(Pseudotime=pData(CDS)$Pseudotime,AUC_cell_data)
#整理后的数据框由于列名格式发生了变化，这里需要重新修改，以便于后续图形的图例可以正常显示regulons名称
colnames(AUC_cell_data1) <- c("Pseudotime",colnames(AUC_cell_data))
#将宽数据整理成长数据，方便使用ggplot2绘制多曲线图
AUC_cell_data2=data.frame()
for (i in 1:(ncol(AUC_cell_data1)-1)){
  a=data.frame(cell = rownames(AUC_cell_data1),Pseudotime=AUC_cell_data1[,1],
               ACU_score=AUC_cell_data1[,1+i],regulons=rep(colnames(AUC_cell_data1)[1+i],nrow(AUC_cell_data1)))
  AUC_cell_data2=rbind(AUC_cell_data2,a)
} 

##7.2 ggplot2绘制图形
##绘制非线性拟合曲线，将每个细胞坐落在拟时间轴上的点拟合成曲线
#按照指定顺序绘制不同regulons曲线
AUC_cell_data2$regulons <- factor(AUC_cell_data2$regulons,levels = unique(AUC_cell_data2$regulons),ordered = T)
#数据映射到图形上
p1 <- ggplot(AUC_cell_data2)+geom_smooth(aes(Pseudotime,ACU_score,color=regulons),method = "loess",se = F,size = 1.2);p1
color <- colorRampPalette(c("#DC143C","#0000FF","#00BFFF","#7FFF00","#FF0000"))(length(unique(AUC_cell_data2$regulons)))
p2 <- p1+scale_color_manual(values = color);p2
p3 <- p2+labs(color="TF Regulons",x="Pseudo-time",y="AUC score",
              title="Pseudo-time curve of TF regulons");p3
mytheme<-theme_bw()+theme(text = element_text(family = "sans"),  #调整图形字体型号
                          plot.title = element_text(size = rel(1.2),hjust = 0.5),  #图形标题字体大小及居中
                          axis.title = element_text(size = rel(1)),  #坐标轴标题字体大小
                          axis.text = element_text(size=rel(0.8),colour = "black"),  #坐标轴刻度的字体大小
                          legend.text = element_text(size = rel(0.8)),  #图例标题字体大小
                          legend.title = element_text(size = rel(1)),  #图例内容字体大小
                          plot.margin=unit(x=c(0.2,0.2,0.2,0.2),units="inches"),  #图形外边框的间距
                          panel.grid = element_blank()) #去除网格线
p4<-p3+mytheme;p4
ggsave(p4,filename = "curve_regulons_all.pdf",width = 18,height = 15,units = "cm")

##7.3 挑选几个regulons绘制拟时间曲线
AUC_cell_data3=subset(AUC_cell_data2,regulons == c("MYC_extended (22g)","SPI1_extended (334g)","CREM (12g)","IRF7 (354g)"))
AUC_cell_data3$regulons <- factor(AUC_cell_data3$regulons,levels = unique(AUC_cell_data3$regulons),ordered = T)
p1 <- ggplot(AUC_cell_data3)+geom_smooth(aes(Pseudotime,ACU_score,color=regulons),method = "loess",se = F,size = 1.2);p1
color <- colorRampPalette(c("#DC143C","#0000FF","#00BFFF"))(length(unique(AUC_cell_data3$regulons)))
p2 <- p1+scale_color_manual(values = color);p2
p3 <- p2+labs(color="TF Regulons",x="Pseudo-time",y="AUC score",
              title="Pseudo-time curve of TF regulons");p3
mytheme<-theme_bw()+theme(text = element_text(family = "sans"),  #调整图形字体型号
                          plot.title = element_text(size = rel(1.2),hjust = 0.5),  #图形标题字体大小及居中
                          axis.title = element_text(size = rel(1)),  #坐标轴标题字体大小
                          axis.text = element_text(size=rel(0.8),colour = "black"),  #坐标轴刻度的字体大小
                          legend.text = element_text(size = rel(0.8)),  #图例标题字体大小
                          legend.title = element_text(size = rel(1)),  #图例内容字体大小
                          plot.margin=unit(x=c(0.2,0.2,0.2,0.2),units="inches"),  #图形外边框的间距
                          panel.grid = element_blank()) #去除网格线
p4<-p3+mytheme;p4
ggsave(p4,filename = "curve_regulons_4.pdf",width = 18,height = 15,units = "cm")

```

[scihub2](https://tool.yovisun.com/scihub/)<br/>
[wiki百科](https://zh.wikipedia.org/wiki/Wikipedia:%E9%A6%96%E9%A1%B5)<br/>

# RNA 数据库条目
### snoRNA/snRNA database
刚开始的时候，snoRNA和snRNA没有明显的差别，可以看到很多早期文献还在称U3为snRNA。所以合并了snRNA/snoRNA。<br/>
[酵母的snoRNA数据库](https://people.biochem.umass.edu/fournierlab/snornadb/mastertable.php)<br/>
[人的snoRNA数据库](https://www-snorna.biotoul.fr/getseq.php)。可以从这个网站下载snoRNA的序列。<br/>


### 人类基因组
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



#### 一些软件
```
python38：
hicexplorer # pip install hicexplorer
hicrep # python版hicrep
```