library(Signac)
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DoubletFinder)
set.seed(1234)
x="MOSCCA0019-all"

RNAs <- Read10X("/data/home/ruanlab/huangxingyu/Haoxi20230215/mouse_brain_markby_process2/MOSCCA0019-all/MOSCCA0019-all_SI-NA-D6/10xGEX/outs/raw_feature_bc_matrix")
ATACs <- Read10X_h5("raw_peak_bc_matrix.h5")
PETs = read.table("/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/4.calPETtestrmblacklist/MOSCCA0019-all.raw.PET.rb.bedpe.gz.cg.stat")
PETflts = PETs[PETs$PET_intragt1000>=50,]
RNA_BC <- colnames(RNAs)
ATAC_BC <- colnames(ATACs)
PET_BC <- rownames(PETflts)
cor_BC <- intersect(intersect(RNA_BC,ATAC_BC),PET_BC)
print(sprintf("There are %d RNA barcodes, there are %d ATAC barcodes, there are %d PET barcodes, and there are %d shared barcodes", length(RNA_BC), length(ATAC_BC), nrow(PETs), length(cor_BC)))
cor_RNAs <- RNAs[,cor_BC]
cor_ATACs <- ATACs[,cor_BC]
cor_PETs <- PETs[cor_BC,]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels
genome(annotations) <- "mm10"

## RNA
mydat <- CreateSeuratObject(cor_RNAs, project = "brain", min.cells = 5) # min.features = 50
mydat[['percent.mt']] <- PercentageFeatureSet(object = mydat, pattern = "^mt-")
mydat[['percent.ribo']] <- PercentageFeatureSet(object = mydat, pattern = "^Rp[sl]")
mydat <- mydat[!grepl("^mt-", rownames(mydat)), ]
mydat <- mydat[!grepl("^Rp[sl]", rownames(mydat)), ]
mydat <- mydat[!grepl("Gm42418", rownames(mydat)), ]

## ATAC
grange.counts <- StringToGRanges(rownames(cor_ATACs), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
cor_ATACs <- cor_ATACs[as.vector(grange.use),]
chrom_assay <- CreateChromatinAssay(counts = cor_ATACs, sep = c(":", "-"), genome = 'mm10', fragments = "/data/home/ruanlab/huangxingyu/Haoxi20230215/mouse_brain_markby_process2/MOSCCA0019-all/MOSCCA0019-all_SI-NA-D6/bulkChiatac/scATACpro/atac_fragments.tsv.gz", min.cells = 5, annotation = annotations)
mydat[["ATAC"]] <- chrom_assay
DefaultAssay(mydat) <- "ATAC"
mydat <- TSSEnrichment(mydat)
mydat <- NucleosomeSignal(mydat)
mydat$blacklist_fraction <- FractionCountsInRegion(object = mydat, assay = 'ATAC', regions = blacklist_mm10)

# filter
mydatflt <- subset(mydat, subset = nCount_RNA > 200 & nCount_ATAC > 200 & percent.mt < 20 & percent.ribo < 50)
cat(sprintf("There are %d cell\n", ncol(mydatflt)))
mydatflt <- AddMetaData(object = mydatflt, metadata = cor_PETs[rownames(cor_PETs) %in% rownames(mydat@meta.data),])

pdf(paste(x,".count.pdf",sep=""), width=15)
p1 = VlnPlot(mydatflt, features = c("nFeature_RNA", "nCount_RNA", "nFeature_ATAC", "nCount_ATAC", "percent.mt", 'percent.ribo'), pt.size = 0, ncol = 3)
print(p1)
p1 = VlnPlot(mydatflt, features = c('TSS.enrichment','TSS.percentile',"PET_intra", 'PET_inter', "PET_intragt1000", 'PET_intrale1000'), pt.size = 0, ncol = 3)
print(p1)
p1 = VlnPlot(mydatflt, features = c("nFeature_RNA", "nCount_RNA", "nFeature_ATAC", "nCount_ATAC", "percent.mt", 'percent.ribo'), pt.size = 0, log = TRUE, ncol = 3)
print(p1)
p1 = VlnPlot(mydatflt, features = c('TSS.enrichment','TSS.percentile',"PET_intra", 'PET_inter', "PET_intragt1000", 'PET_intrale1000'), pt.size = 0, log = TRUE, ncol = 3)
print(p1)
p1 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "percent.ribo")
p3 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1 + p2 + p3)
p1 = FeatureScatter(mydatflt, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC")
p2 = FeatureScatter(mydatflt, feature1 = "nCount_ATAC", feature2 = "nCount_RNA")
p3 = FeatureScatter(mydatflt, feature1 = "nCount_ATAC", feature2 = "nFeature_RNA")
print(p1 + p2 + p3)
p1 = FeatureScatter(mydatflt, feature1 = "PET_intra", feature2 = "PET_inter")
p2 = FeatureScatter(mydatflt, feature1 = "PET_intra", feature2 = "PET_intrale1000")
p3 = FeatureScatter(mydatflt, feature1 = "PET_intra", feature2 = "nCount_ATAC")
print(p1 + p2 + p3)
dev.off()

mydatflt <- subset(mydat, subset = nCount_RNA > 400 & nCount_ATAC > 2000 & percent.mt < 20 & percent.ribo < 50)
cat(sprintf("There are %d cell\n", ncol(mydatflt)))
mydatflt <- AddMetaData(object = mydatflt, metadata = cor_PETs[rownames(cor_PETs) %in% rownames(mydat@meta.data),])
saveRDS(mydatflt, file="MOSCCA0019-all.rds")

pdf(paste(x,".test1.pdf",sep=""))
DensityScatter(mydatflt, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE) + ylim(c(0,10))
dev.off()

# ATAC analysis
mydatflt <- RunTFIDF(mydatflt)
mydatflt <- FindTopFeatures(mydatflt, min.cutoff = 'q0')
mydatflt <- RunSVD(mydatflt)
mydatflt <- RunUMAP(object = mydatflt, reduction = 'lsi', dims = 2:30, reduction.name = "lsiumap", reduction.key = "lsiUMAP_")
mydatflt <- RunTSNE(object = mydatflt, reduction = 'lsi', dims = 2:30, reduction.name = "lsitsne", reduction.key = "lsiSNE_")
mydatflt <- FindNeighbors(object = mydatflt, reduction = 'lsi', dims = 2:30)
mydatflt <- FindClusters(object = mydatflt, verbose = FALSE, algorithm = 3)

# RNA analysis
DefaultAssay(mydatflt) = "RNA"
mydatflt <- NormalizeData(object = mydatflt)
mydatflt <- FindVariableFeatures(object = mydatflt,nfeatures = 3000)
mydatflt <- ScaleData(object = mydatflt, features = rownames(mydatflt))
mydatflt <- RunPCA(mydatflt, features = VariableFeatures(mydatflt))
mydatflt <- FindNeighbors(object = mydatflt, dims = 1:30)
mydatflt <- FindClusters(mydatflt, resolution = 0.8)
mydatflt <- RunTSNE(mydatflt, dims = 1:30, reduction.name="pcatsne", reduction.key = "pcaTSNE_", check_duplicates = FALSE) 
mydatflt <- RunUMAP(mydatflt, dims = 1:30, reduction.name="pcaumap", reduction.key = "pcaUMAP_")

# WNN analysis
mydatflt <- FindMultiModalNeighbors(mydatflt, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
mydatflt <- FindClusters(mydatflt, graph.name = "wsnn", algorithm = 3, resolution = 0.5)
mydatflt <- RunTSNE(mydatflt, nn.name = "weighted.nn", reduction.name = "wnntsne", reduction.key = "wnnTSNE_")
mydatflt <- RunUMAP(mydatflt, nn.name = "weighted.nn", reduction.name = "wnnumap", reduction.key = "wnnUMAP_")

pdf(paste(x,".test2.pdf",sep=""),width=15)
p1 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsitsne", group.by="ATAC_snn_res.0.8") + ggtitle("lsitsne") + NoLegend()
p2 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsiumap", group.by="ATAC_snn_res.0.8") + ggtitle("lsiumap") + NoLegend()

p3 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcatsne", group.by="RNA_snn_res.0.8") + ggtitle("pcatsne") + NoLegend()
p4 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcaumap", group.by="RNA_snn_res.0.8") + ggtitle("pcaumap") + NoLegend()

p5 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "wnntsne", group.by="wsnn_res.0.5") + ggtitle("wnntsne") + NoLegend()
p6 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "wnnumap", group.by="wsnn_res.0.5") + ggtitle("wnnumap") + NoLegend()
print((p1 + p3 + p5)/(p2 + p4 + p6))
dev.off()

pdf(paste(x,".test3.pdf",sep=""),width=15,height=20)
FeaturePlot(object = mydatflt, reduction = "lsitsne", features = c('Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Pdgfra', 'Mbp','Apoe', 'Cldn5', 'C1qa', 'Igfbpl1'), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
FeaturePlot(object = mydatflt, reduction = "lsiumap", features = c('Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Pdgfra', 'Mbp','Apoe', 'Cldn5', 'C1qa', 'Igfbpl1'), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
FeaturePlot(object = mydatflt, reduction = "pcatsne", features = c('Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Pdgfra', 'Mbp','Apoe', 'Cldn5', 'C1qa', 'Igfbpl1'), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
FeaturePlot(object = mydatflt, reduction = "pcaumap", features = c('Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Pdgfra', 'Mbp','Apoe', 'Cldn5', 'C1qa', 'Igfbpl1'), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
FeaturePlot(object = mydatflt, reduction = "wnntsne", features = c('Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Pdgfra', 'Mbp','Apoe', 'Cldn5', 'C1qa', 'Igfbpl1'), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
FeaturePlot(object = mydatflt, reduction = "wnnumap", features = c('Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Pdgfra', 'Mbp','Apoe', 'Cldn5', 'C1qa', 'Igfbpl1'), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
dev.off()
saveRDS(mydatflt, file="MOSCCA0019-all.final.rds")

pdf(paste(x,".test4.pdf",sep=""),width=15)
p1 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsitsne", group.by="RNA_snn_res.0.8") + ggtitle("lsitsne") + NoLegend()
p2 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsiumap", group.by="RNA_snn_res.0.8") + ggtitle("lsiumap") + NoLegend()

p3 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcatsne", group.by="RNA_snn_res.0.8") + ggtitle("pcatsne") + NoLegend()
p4 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcaumap", group.by="RNA_snn_res.0.8") + ggtitle("pcaumap") + NoLegend()

p5 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "wnntsne", group.by="RNA_snn_res.0.8") + ggtitle("wnntsne") + NoLegend()
p6 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "wnnumap", group.by="RNA_snn_res.0.8") + ggtitle("wnnumap") + NoLegend()
print((p1 + p3 + p5)/(p2 + p4 + p6))
dev.off()

write.table(mydatflt@meta.data, file=paste0(x,".metadata.txt"),sep="\t",quote=F)


# 寻找最优pK值
pcselect = 30
sweep = paramSweep_v3(mydatflt, PCs=1:pcselect, sct=F)
sweep.stats <- summarizeSweep(sweep, GT = FALSE)
pdf("pk.test.pdf")
bcmvn <- find.pK(sweep.stats)
dev.off()
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()#提取最佳pk值

# 检测双细胞
# DoubletRate = 0.054 # 直接查表，7000细胞对应的doubletsrate是~5.4%
# DoubletRate = ncol(mydatflt)*8*1e-6 # 按每增加1808个细胞，双细胞比率增加千分之8来计算
DoubletRate = ncol(mydatflt)*8*1e-6 # 更通用#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat clusters中来混双细胞
homotypic.prop <- modelHomotypic(mydatflt$RNA_snn_res.0.8) # 最好提供celltype, 而不是seurat cluster#计算双细胞比例
nExp_poi <- round(DoubletRate * ncol(mydatflt)) # 使用同源双细胞比例对计算的双细胞比例进行校正
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))
##使用确定好的参数鉴定doublets
mydatflt <- doubletFinder_v3(mydatflt, PCs = 1:pcselect, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN=F, sct =F)
## 结果展示，分类结果在mydatflt@meta.data中
DF.name = colnames(mydatflt@meta.data)[grepl("DF.classification", colnames(mydatflt@meta.data))]
pANN.name = colnames(mydatflt@meta.data)[grepl("pANN", colnames(mydatflt@meta.data))]
mydatflt <- AddMetaData(mydatflt, mydatflt[[pANN.name]], col.name = "RNApANN")
mydatflt <- AddMetaData(mydatflt, mydatflt[[DF.name]], col.name = "RNADF.classifications")
mydatflt[[pANN.name]] <- NULL
mydatflt[[DF.name]] <- NULL

pdf(paste(x,".test5.doublet.pdf",sep=""), width=15)
p1 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsitsne", group.by = c("RNADF.classifications", "RNA_snn_res.0.8")) + ggtitle("DFRNAlsitsne") + NoLegend()
p2 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsiumap", group.by = c("RNADF.classifications", "RNA_snn_res.0.8")) + ggtitle("DFRNAlsiumap") + NoLegend()

p3 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcatsne", group.by = c("RNADF.classifications", "RNA_snn_res.0.8")) + ggtitle("DFRNApcatsne") + NoLegend()
p4 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcaumap", group.by = c("RNADF.classifications", "RNA_snn_res.0.8")) + ggtitle("DFRNApcaumap") + NoLegend()

p5 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "wnntsne", group.by = c("RNADF.classifications", "RNA_snn_res.0.8")) + ggtitle("DFRNAwnntsne") + NoLegend()
p6 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "wnnumap", group.by = c("RNADF.classifications", "RNA_snn_res.0.8")) + ggtitle("DFRNAwnnumap") + NoLegend()
print((p1 + p3 + p5)/(p2 + p4 + p6))
dev.off()

saveRDS(mydatflt, file="MOSCCA0019-all.doublet.rds")
write.table(mydatflt@meta.data, file=paste0(x,".doublet.metadata.txt"),sep="\t",quote=F)

mydatflt <- subset(mydatflt, subset = RNADF.classifications == "Singlet")
mydatflt[["wsnn_res.0.5"]] <- NULL
mydatflt[["ATAC.weight"]] <- NULL
mydatflt[["RNA.weight"]] <- NULL
mydatflt[["RNA_snn_res.0.8"]] <- NULL
mydatflt[["seurat_clusters"]] <- NULL
mydatflt[["ATAC_snn_res.0.8"]] <- NULL
mydatflt@reductions$lsi = NULL
mydatflt@reductions$lsiumap = NULL
mydatflt@reductions$lsitsne = NULL
mydatflt@reductions$pca = NULL
mydatflt@reductions$pcaumap = NULL
mydatflt@reductions$pcatsne = NULL
mydatflt@reductions$wnnumap = NULL
mydatflt@reductions$wnntsne = NULL

# RNA analysis
DefaultAssay(mydatflt) = "RNA"
mydatflt <- NormalizeData(object = mydatflt)
mydatflt <- FindVariableFeatures(object = mydatflt,nfeatures = 3000)
mydatflt <- ScaleData(object = mydatflt, features = rownames(mydatflt))
mydatflt <- RunPCA(mydatflt, features = VariableFeatures(mydatflt))
mydatflt <- FindNeighbors(object = mydatflt, dims = 1:30)
mydatflt <- FindClusters(mydatflt, resolution = 0.8)
mydatflt <- RunTSNE(mydatflt, dims = 1:30, reduction.name="pcatsne", reduction.key = "pcaTSNE_", check_duplicates = FALSE) 
mydatflt <- RunUMAP(mydatflt, dims = 1:30, reduction.name="pcaumap", reduction.key = "pcaUMAP_")

# ATAC analysis
DefaultAssay(mydatflt) = "ATAC"
mydatflt <- RunTFIDF(mydatflt)
mydatflt <- FindTopFeatures(mydatflt, min.cutoff = 'q0')
mydatflt <- RunSVD(mydatflt)
mydatflt <- RunUMAP(object = mydatflt, reduction = 'lsi', dims = 2:30, reduction.name = "lsiumap", reduction.key = "lsiUMAP_")
mydatflt <- RunTSNE(object = mydatflt, reduction = 'lsi', dims = 2:30, reduction.name = "lsitsne", reduction.key = "lsiSNE_")
mydatflt <- FindNeighbors(object = mydatflt, reduction = 'lsi', dims = 2:30)
mydatflt <- FindClusters(object = mydatflt, verbose = FALSE, algorithm = 3)

# PET analysis
pet = readRDS("/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/4.calPETtestrmblacklistflt/MOSCCA0019-all.1000k.bedpe.mat.original.PET.rds")
mydatflt[["PET"]] <- CreateAssayObject(pet)
DefaultAssay(mydatflt) <- "PET"
mydatflt <- FindTopFeatures(mydatflt, assay="PET", min.cutoff = "q1")
mydatflt = RunSVD(mydatflt, assay = "PET", reduction.key = "PETLSI_", reduction.name = "petlsi")
pdf(paste(x,".test5.DepthCor.pdf",sep=""), width=15)
DepthCor(mydatflt, assay = "PET", reduction = "petlsi")
dev.off()
mydatflt <- FindNeighbors(mydatflt, reduction = 'petlsi', dims = 1:20, k.param = 30)
mydatflt <- FindClusters(mydatflt, resolution = 0.5, algorithm = 3)
mydatflt <- RunUMAP(mydatflt, reduction = "petlsi", dims = 1:20, reduction.name = "petlsiumap", reduction.key = "petlsiUMAP_")
mydatflt <- RunTSNE(mydatflt, reduction = "petlsi", dims = 1:20, reduction.name = "petlsitsne", reduction.key = "petlsiTSNE_", check_duplicates = FALSE)

# wnn 
mydatflt <- FindMultiModalNeighbors(mydatflt, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30), snn.graph.name = "pcalsi", weighted.nn.name = "weightednn_pcalsi")
mydatflt <- RunUMAP(mydatflt, nn.name = "weightednn_pcalsi", reduction.name = "pcalsiumap", reduction.key = "pcalsiUMAP_")

mydatflt <- FindMultiModalNeighbors(mydatflt, reduction.list = list("pca", "petlsi"), dims.list = list(1:30, 1:20), snn.graph.name = "pcapetlsi", weighted.nn.name = "weightednn_pcapetlsi")
mydatflt <- RunUMAP(mydatflt, nn.name = "weightednn_pcapetlsi", reduction.name = "pcapetlsiump", reduction.key = "pcapetlsiUMAP_")

mydatflt <- FindMultiModalNeighbors(mydatflt, reduction.list = list("lsi", "petlsi"), dims.list = list(2:30, 1:20), snn.graph.name = "lsipetlsi", weighted.nn.name = "weightednn_lsipetlsi")
mydatflt <- RunUMAP(mydatflt, nn.name = "weightednn_lsipetlsi", reduction.name = "lsipetlsiumap", reduction.key = "lsipetlsiUMAP_")

mydatflt <- FindMultiModalNeighbors(mydatflt, reduction.list = list("pca", "lsi", "petlsi"), dims.list = list(1:30, 2:30, 1:20), snn.graph.name = "wsnn3", weighted.nn.name = "weightednn3")
mydatflt <- RunUMAP(mydatflt, nn.name = "weightednn3", reduction.name = "weightednn3umap", reduction.key = "weightednn3UMAP_")

pdf(paste(x,".wnn3all.pdf",sep=""),width=15)
p1 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsitsne", group.by="RNA_snn_res.0.8") + ggtitle("lsitsne") + NoLegend()
p2 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsiumap", group.by="RNA_snn_res.0.8") + ggtitle("lsiumap") + NoLegend()
p3 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcatsne", group.by="RNA_snn_res.0.8") + ggtitle("pcatsne") + NoLegend()
p4 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcaumap", group.by="RNA_snn_res.0.8") + ggtitle("pcaumap") + NoLegend()
p5 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "petlsitsne", group.by="RNA_snn_res.0.8") + ggtitle("petlsitsne") 
p6 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "petlsiumap", group.by="RNA_snn_res.0.8") + ggtitle("petlsiumap") 
print(p1 + p3 + p5)
print(p2 + p4 + p6)

p2 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcalsiumap", group.by="RNA_snn_res.0.8") + ggtitle("pcalsiumap") + NoLegend()
p4 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "pcapetlsiump", group.by="RNA_snn_res.0.8") + ggtitle("pcapetlsiump") + NoLegend()
p6 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "lsipetlsiumap", group.by="RNA_snn_res.0.8") + ggtitle("lsipetlsiumap") 
p8 = DimPlot(object = mydatflt, pt.size = 0.1, reduction = "weightednn3umap", group.by="RNA_snn_res.0.8") + ggtitle("weightednn3umap") + NoLegend()
print(p2 + p4 + p6)
print(p8)
dev.off()
saveRDS(mydatflt, file="MOSCCA0019-all.doublet.cluster.wnn.rds")

# annotation
rnadata = mydatflt
DefaultAssay(rnadata) = "RNA"
rnadata[["ATAC"]] = NULL
rnadata[["PET"]] = NULL
rnadata@reductions$lsi = NULL
rnadata@reductions$lsiumap = NULL
rnadata@reductions$lsitsne = NULL
rnadata@reductions$petlsi = NULL
rnadata@reductions$petlsiumap = NULL
rnadata@reductions$petlsitsne = NULL

#### cca
reference=readRDS("/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/0.referenceMouseBrainMap/mouseBrain.brainPart.reference.rds")
cca.peak.anchor <- FindIntegrationAnchors(object.list = c(rnadata, reference), anchor.features = 3000, reduction = "cca")
saveRDS(cca.peak.anchor, "MOSCCA0019-all.cca.peak.anchor.rds")
cca.Integrate <- IntegrateData(anchorset = cca.peak.anchor)
cca.Integrate <- ScaleData(cca.Integrate, features=rownames(cca.Integrate))
cca.Integrate <- RunPCA(cca.Integrate, features = VariableFeatures(cca.Integrate))
cca.Integrate <- RunTSNE(cca.Integrate, dims = 1:30, reduction="pca")
cca.Integrate <- RunUMAP(cca.Integrate, dims = 1:30, reduction="pca")

pdf("MOSCCA0019-all.cca.peak.plot.pdf", height=10, width=20)
p1 = DimPlot(cca.Integrate, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group.by = "orig.ident") + ggtitle("pca")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group.by = "orig.ident") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group.by = "orig.ident") + ggtitle("umap")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, split.by = "orig.ident") + ggtitle("pca")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, split.by = "orig.ident") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, split.by = "orig.ident") + ggtitle("umap")
print(p1)
dev.off()
saveRDS(cca.Integrate, file = "MOSCCA0019-all.cca.Integrate.all.rds")


##### harmony
library(harmony)
harmonymerge <- merge(rnadata, y = c(reference), project = "mouseBrain")
harmonymerge <- NormalizeData(object = harmonymerge)
harmonymerge <- FindVariableFeatures(object = harmonymerge, nfeatures = 3000)
harmonymerge <- ScaleData(harmonymerge, features=rownames(harmonymerge))
harmonymerge <- RunPCA(harmonymerge, features = VariableFeatures(object = harmonymerge))

harmonymerge <- RunHarmony(harmonymerge, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "rnaharmonypca", max.iter.harmony = 100)
harmonymerge@reductions$pca = harmonymerge@reductions$rnaharmonypca
harmonymerge <- RunUMAP(harmonymerge, dims = 1:30, reduction="rnaharmonypca", reduction.name="umapharmony", reduction.key = "tsneharmonyUMAP_")
pdf("MOSCCA0019-all.harmony1.pdf", height=10, width=20)
p1 = DimPlot(harmonymerge, reduction ="umapharmony", label = TRUE, label.size=2, raster=FALSE, group.by = "orig.ident") + ggtitle("umapharmony")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="umapharmony", raster=FALSE, group.by = "orig.ident", split.by = "orig.ident") + ggtitle("umapharmony")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="umapharmony", label = TRUE, raster=FALSE, group.by = "celltype", split.by = "orig.ident") + ggtitle("umapharmony")
print(p1)
dev.off()
saveRDS(harmonymerge, file = "MOSCCA0019-all.harmony.all.rds")


##### transfer
cca.transfer.anchors <- FindTransferAnchors(reference = reference, query = rnadata, dims = 1:30, reduction = 'cca' )
saveRDS(cca.transfer.anchors, "MOSCCA0019-all.cca.transfer.anchors.rds")
cca.predicted.labels <- TransferData(anchorset = cca.transfer.anchors, refdata = reference$celltype, weight.reduction = rnadata[['pca']], dims = 1:30)
rnadata <- AddMetaData(object = rnadata, metadata = cca.predicted.labels)
rnadata <- AddMetaData(object = rnadata, metadata = rnadata$predicted.id, col.name = 'celltype')
saveRDS(rnadata, file="MOSCCA0019-all.cca.transfer.rds")
for (x in grep("RNA_snn_res.",colnames(SCG@meta.data),value = T)){SCG[[x]] <- NULL}



rnadata = mydatflt
DefaultAssay(rnadata) = "RNA"
rnadata[["ATAC"]] = NULL
rnadata[["PET"]] = NULL
rnadata@reductions$lsi = NULL
rnadata@reductions$lsiumap = NULL
rnadata@reductions$lsitsne = NULL
rnadata@reductions$petlsi = NULL
rnadata@reductions$petlsiumap = NULL
rnadata@reductions$petlsitsne = NULL

#### cca
reference=readRDS("/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/8.mergeSampleBrainFinal/mouseBrain.P365.rds")
reference@assays$ATAC=NULL
reference@assays$PET=NULL

cca.peak.anchor <- FindIntegrationAnchors(object.list = c(rnadata, reference), anchor.features = 3000, reduction = "cca")
saveRDS(cca.peak.anchor, "MOSCCA0019-all.P365.cca.peak.anchor.rds")
cca.Integrate <- IntegrateData(anchorset = cca.peak.anchor)
cca.Integrate <- ScaleData(cca.Integrate, features=rownames(cca.Integrate))
cca.Integrate <- RunPCA(cca.Integrate, features = VariableFeatures(cca.Integrate))
cca.Integrate <- RunUMAP(cca.Integrate, dims = 1:30, reduction="pca")
cca.Integrate <- RunTSNE(cca.Integrate, dims = 1:30, reduction="pca")
cca.Integrate@meta.data[is.na(cca.Integrate$Age),]$Age = "newP365"

pdf("MOSCCA0019-all.P365.cca.peak.plot.pdf", height=10, width=20)
p1 = DimPlot(cca.Integrate, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group.by = "Age") + ggtitle("pca")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group.by = "Age") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group.by = "Age") + ggtitle("umap")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, split.by = "Age") + ggtitle("pca")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, split.by = "Age") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, split.by = "Age") + ggtitle("umap")
print(p1)
dev.off()
saveRDS(cca.Integrate, file = "MOSCCA0019-all.P365.cca.Integrate.all.rds")


harmonymerge <- merge(rnadata, y = c(reference), project = "mouseBrain")
harmonymerge <- NormalizeData(object = harmonymerge)
harmonymerge <- FindVariableFeatures(object = harmonymerge, nfeatures = 3000)
harmonymerge <- ScaleData(harmonymerge, features=rownames(harmonymerge))
harmonymerge <- RunPCA(harmonymerge, features = VariableFeatures(object = harmonymerge))

harmonymerge <- RunHarmony(harmonymerge, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "rnaharmonypca", max.iter.harmony = 100)
harmonymerge@reductions$pca = harmonymerge@reductions$rnaharmonypca
harmonymerge <- RunUMAP(harmonymerge, dims = 1:30, reduction="rnaharmonypca", reduction.name="umapharmony", reduction.key = "tsneharmonyUMAP_")
harmonymerge@meta.data[is.na(harmonymerge$Age),]$Age = "newP365"

pdf("MOSCCA0019-all.harmony2.pdf", height=10, width=20)
p1 = DimPlot(harmonymerge, reduction ="umapharmony", label = TRUE, label.size=2, raster=FALSE, group.by = "Age") + ggtitle("umapharmony")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="umapharmony", raster=FALSE, group.by = "Age", split.by = "Age") + ggtitle("umapharmony")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="umapharmony", label = TRUE, raster=FALSE, group.by = "celltype", split.by = "Age") + ggtitle("umapharmony")
print(p1)
dev.off()
saveRDS(harmonymerge, file = "MOSCCA0019-all.P365.harmony.all.rds")


cca.transfer.anchors <- FindTransferAnchors(reference = reference, query = rnadata, dims = 1:30, reduction = 'cca' )
saveRDS(cca.transfer.anchors, "MOSCCA0019-all.P365.cca.transfer.anchors.rds")
cca.predicted.labels <- TransferData(anchorset = cca.transfer.anchors, refdata = reference$celltype, weight.reduction = rnadata[['pca']], dims = 1:30)
rnadata <- AddMetaData(object = rnadata, metadata = cca.predicted.labels)
rnadata <- AddMetaData(object = rnadata, metadata = rnadata$predicted.id, col.name = 'celltype')
saveRDS(rnadata, file="MOSCCA0019-all.P365.cca.transfer.rds")



map = read.table("/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/celltype.txt",header=T,sep="\t")
mapflt = map[map$group!="combine",]
mydat = readRDS(file="MOSCCA0019-all.doublet.cluster.wnn.rds")
rnadata1 = readRDS(file="MOSCCA0019-all.cca.transfer.rds")
rnadata2 = readRDS(file="MOSCCA0019-all.P365.cca.transfer.rds")
result1 <- left_join(rnadata1@meta.data, mapflt, by = "celltype")
result2 <- left_join(rnadata2@meta.data, mapflt, by = "celltype")

mydat$celltype1 = rnadata1$celltype
mydat$celltype2 = rnadata2$celltype
mydat <- AddMetaData(object = mydat, metadata = result1$bigClass1, col.name = 'bigClass1')
mydat <- AddMetaData(object = mydat, metadata = result2$bigClass1, col.name = 'bigClass2')

pdf(paste(x,".wnn3all.celltype.pdf",sep=""),width=20)
p1 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsitsne", group.by="celltype1") + ggtitle("lsitsne") + NoLegend()
p2 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsiumap", group.by="celltype1") + ggtitle("lsiumap") + NoLegend()
p3 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcatsne", group.by="celltype1") + ggtitle("pcatsne") + NoLegend()
p4 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcaumap", group.by="celltype1") + ggtitle("pcaumap") + NoLegend()
p5 = DimPlot(object = mydat, pt.size = 0.1, reduction = "petlsitsne", group.by="celltype1", label.size = 2) + ggtitle("petlsitsne") 
p6 = DimPlot(object = mydat, pt.size = 0.1, reduction = "petlsiumap", group.by="celltype1", label.size = 2) + ggtitle("petlsiumap") 
print(p1 + p3 + p5)
print(p2 + p4 + p6)
p2 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcalsiumap", group.by="celltype1") + ggtitle("pcalsiumap") + NoLegend()
p4 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcapetlsiump", group.by="celltype1") + ggtitle("pcapetlsiump") + NoLegend()
p6 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsipetlsiumap", group.by="celltype1", label.size = 2) + ggtitle("lsipetlsiumap") 
p8 = DimPlot(object = mydat, pt.size = 0.1, reduction = "weightednn3umap", group.by="celltype1") + ggtitle("weightednn3umap")
print(p2 + p4 + p6)
print(p8)

p1 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsitsne", group.by="celltype2") + ggtitle("lsitsne") + NoLegend()
p2 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsiumap", group.by="celltype2") + ggtitle("lsiumap") + NoLegend()
p3 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcatsne", group.by="celltype2") + ggtitle("pcatsne") + NoLegend()
p4 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcaumap", group.by="celltype2") + ggtitle("pcaumap") + NoLegend()
p5 = DimPlot(object = mydat, pt.size = 0.1, reduction = "petlsitsne", group.by="celltype2", label.size = 2) + ggtitle("petlsitsne") 
p6 = DimPlot(object = mydat, pt.size = 0.1, reduction = "petlsiumap", group.by="celltype2", label.size = 2) + ggtitle("petlsiumap") 
print(p1 + p3 + p5)
print(p2 + p4 + p6)
p2 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcalsiumap", group.by="celltype2") + ggtitle("pcalsiumap") + NoLegend()
p4 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcapetlsiump", group.by="celltype2") + ggtitle("pcapetlsiump") + NoLegend()
p6 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsipetlsiumap", group.by="celltype2", label.size = 2) + ggtitle("lsipetlsiumap") 
p8 = DimPlot(object = mydat, pt.size = 0.1, reduction = "weightednn3umap", group.by="celltype2", label.size = 2) + ggtitle("weightednn3umap")
print(p2 + p4 + p6)
print(p8)
dev.off()
pdf(paste(x,".wnn3all.bigClass.pdf",sep=""),width=20)
p1 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsitsne", group.by="bigClass1") + ggtitle("lsitsne") + NoLegend()
p2 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsiumap", group.by="bigClass1") + ggtitle("lsiumap") + NoLegend()
p3 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcatsne", group.by="bigClass1") + ggtitle("pcatsne") + NoLegend()
p4 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcaumap", group.by="bigClass1") + ggtitle("pcaumap") + NoLegend()
p5 = DimPlot(object = mydat, pt.size = 0.1, reduction = "petlsitsne", group.by="bigClass1", label.size = 2) + ggtitle("petlsitsne") 
p6 = DimPlot(object = mydat, pt.size = 0.1, reduction = "petlsiumap", group.by="bigClass1", label.size = 2) + ggtitle("petlsiumap") 
print(p1 + p3 + p5)
print(p2 + p4 + p6)
p2 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcalsiumap", group.by="bigClass1") + ggtitle("pcalsiumap") + NoLegend()
p4 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcapetlsiump", group.by="bigClass1") + ggtitle("pcapetlsiump") + NoLegend()
p6 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsipetlsiumap", group.by="bigClass1", label.size = 2) + ggtitle("lsipetlsiumap") 
p8 = DimPlot(object = mydat, pt.size = 0.1, reduction = "weightednn3umap", group.by="bigClass1") + ggtitle("weightednn3umap")
print(p2 + p4 + p6)
print(p8)

p1 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsitsne", group.by="bigClass2") + ggtitle("lsitsne") + NoLegend()
p2 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsiumap", group.by="bigClass2") + ggtitle("lsiumap") + NoLegend()
p3 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcatsne", group.by="bigClass2") + ggtitle("pcatsne") + NoLegend()
p4 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcaumap", group.by="bigClass2") + ggtitle("pcaumap") + NoLegend()
p5 = DimPlot(object = mydat, pt.size = 0.1, reduction = "petlsitsne", group.by="bigClass2", label.size = 2) + ggtitle("petlsitsne") 
p6 = DimPlot(object = mydat, pt.size = 0.1, reduction = "petlsiumap", group.by="bigClass2", label.size = 2) + ggtitle("petlsiumap") 
print(p1 + p3 + p5)
print(p2 + p4 + p6)
p2 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcalsiumap", group.by="bigClass2") + ggtitle("pcalsiumap") + NoLegend()
p4 = DimPlot(object = mydat, pt.size = 0.1, reduction = "pcapetlsiump", group.by="bigClass2") + ggtitle("pcapetlsiump") + NoLegend()
p6 = DimPlot(object = mydat, pt.size = 0.1, reduction = "lsipetlsiumap", group.by="bigClass2", label.size = 2) + ggtitle("lsipetlsiumap") 
p8 = DimPlot(object = mydat, pt.size = 0.1, reduction = "weightednn3umap", group.by="bigClass2", label.size = 2) + ggtitle("weightednn3umap")
print(p2 + p4 + p6)
print(p8)
dev.off()