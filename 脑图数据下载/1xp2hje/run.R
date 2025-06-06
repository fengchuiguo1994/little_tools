library(Seurat)
library(dplyr)
library(ggplot2)
set.seed(1234)

counts <- Read10X(data.dir = "1xp2hje")
scRNA <- CreateSeuratObject(counts = counts, project = "1xp2hje", assay = "RNA", min.cells = 3, min.features = 3)
saveRDS(scRNA, file = "1xp2hje.raw.data.rds")

sce <- Read10X_h5(filename = "umi_counts.h5")
sce <- CreateSeuratObject(counts = sce)
saveRDS(sce, file = "1xp2hje.raw.data.h5.rds")

# length(rownames(scRNA@assays$RNA@counts))
# length(rownames(sce@assays$RNA@counts))
# length(intersect(rownames(sce@assays$RNA@counts),rownames(scRNA@assays$RNA@counts)))

metadata = read.csv("sample_metadata.csv",header=T)
rownames(metadata) = metadata$X
# table(rownames(scRNA@meta.data) == rownames(metadata))
scRNA = AddMetaData(scRNA, metadata)

anno1 = read.csv("cluster.membership.csv", header=T)
names(anno1)[2] = "cluster_id"
anno2 = read.csv("cluster.annotation.csv", header=T)
anno = left_join(anno1,anno2,by="cluster_id")

fltlist = intersect(metadata$X, anno$X)
droplist = setdiff(metadata$X, anno$X)

scRNAflt = subset(scRNA, subset=(X %in% fltlist))
scRNAdrop = subset(scRNA, subset=(X %in% droplist))
# table(rownames(scRNAflt@meta.data) == anno$X)
annod = anno[,-1]
rownames(annod) = anno$X
scRNAflt = AddMetaData(scRNAflt, annod)


scRNAflt <- NormalizeData(object = scRNAflt)
scRNAflt <- FindVariableFeatures(object = scRNAflt, nfeatures = 3000)
scRNAflt <- ScaleData(scRNAflt, features=rownames(scRNAflt))
scRNAflt <- RunPCA(scRNAflt, features = VariableFeatures(object = scRNAflt))
scRNAflt <- RunUMAP(scRNAflt, dims = 1:30)
scRNAflt <- RunTSNE(scRNAflt, dims = 1:30)

pdf("1xp2hje.class.pdf", height=10, width=20)
p1 = DimPlot(scRNAflt, reduction = "umap", raster = FALSE, label = T, group.by = "cluster_label", label.size = 6) + ggtitle("umap")
print(p1)
p1 = DimPlot(scRNAflt, reduction = "umap", raster = FALSE, label = T, group.by = "subclass_label", label.size = 6) + ggtitle("umap")
print(p1)
p1 = DimPlot(scRNAflt, reduction = "umap", raster = FALSE, label = T, group.by = "class_label", label.size = 6) + ggtitle("umap")
print(p1)

p1 = DimPlot(scRNAflt, reduction = "tsne", raster = FALSE, label = T, group.by = "cluster_label", label.size = 6) + ggtitle("tsne")
print(p1)
p1 = DimPlot(scRNAflt, reduction = "tsne", raster = FALSE, label = T, group.by = "subclass_label", label.size = 6) + ggtitle("tsne")
print(p1)
p1 = DimPlot(scRNAflt, reduction = "tsne", raster = FALSE, label = T, group.by = "class_label", label.size = 6) + ggtitle("tsne")
print(p1)
dev.off()

saveRDS(scRNAflt, file = "1xp2hje.cluster.rds")




# aaa = scRNAflt@meta.data
# head(aaa[aaa$subclass_label=="Doublet" | aaa$subclass_label=="Low Quality" | aaa$class_label=="Low Quality",])
# table(aaa[aaa$subclass_label=="Doublet",]$cluster_label)
# table(aaa[aaa$subclass_label=="Low Quality",]$cluster_label)
# table(aaa[aaa$class_label=="Low Quality",]$cluster_label)
SCG = subset(scRNAflt, subset=(subclass_label!="Doublet" & subclass_label!="Low Quality" & class_label!="Low Quality"))
SCG@reductions$pca = NULL
SCG@reductions$rawpca = NULL
SCG@reductions$umap = NULL
SCG@reductions$tsne = NULL

SCG <- NormalizeData(object = SCG)
SCG <- FindVariableFeatures(object = SCG, nfeatures = 3000)
SCG <- ScaleData(SCG, features=rownames(SCG))
SCG <- RunPCA(SCG, features = VariableFeatures(object = SCG))
SCG <- RunUMAP(SCG, dims = 1:30)
SCG <- RunTSNE(SCG, dims = 1:30)

pdf("1xp2hje.final.class.pdf", height=10, width=20)
p1 = DimPlot(SCG, reduction = "umap", raster = FALSE, label = T, group.by = "cluster_label", label.size = 6) + ggtitle("umap")
print(p1)
p1 = DimPlot(SCG, reduction = "umap", raster = FALSE, label = T, group.by = "subclass_label", label.size = 6) + ggtitle("umap")
print(p1)
p1 = DimPlot(SCG, reduction = "umap", raster = FALSE, label = T, group.by = "class_label", label.size = 6) + ggtitle("umap")
print(p1)

p1 = DimPlot(SCG, reduction = "tsne", raster = FALSE, label = T, group.by = "cluster_label", label.size = 6) + ggtitle("tsne")
print(p1)
p1 = DimPlot(SCG, reduction = "tsne", raster = FALSE, label = T, group.by = "subclass_label", label.size = 6) + ggtitle("tsne")
print(p1)
p1 = DimPlot(SCG, reduction = "tsne", raster = FALSE, label = T, group.by = "class_label", label.size = 6) + ggtitle("tsne")
print(p1)
dev.off()
saveRDS(SCG, file = "1xp2hje.final.cluster.rds")

P95 = readRDS("/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/8.mergeSampleBrainFinal/mouseBrain.P95.rds")
aa = P95@meta.data
non = aa[aa$bigClass1!="Excitatory_Neurons" & aa$bigClass1!="Inhibitory_Neurons",]
neuron = aa[aa$bigClass1=="Excitatory_Neurons" | aa$bigClass1=="Inhibitory_Neurons",]
corneu = neuron[neuron$celltype2 == "TEGLU" | neuron$celltype2 == "TEINH" | neuron$celltype2 == "MSN",]
P95$cellid = rownames(aa)
P95flt = subset(P95, subset=(cellid %in% c(rownames(non), rownames(corneu))))
P95 = P95flt
mydat = as.data.frame(table(P95@meta.data$celltype))
flt = as.character(mydat[mydat$Freq >= 50,]$Var1)
P95 = subset(P95, subset = (celltype %in% flt))
P95 = subset(P95, subset = (celltype != "VSMCA"))

P95 = DietSeurat(P95, counts = TRUE, scale.data = F, assays = "RNA")
for (i in c("percent.mt", "percent.ribo", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "TSS.percentile", "nucleosome_signal", "nucleosome_percentile", "blacklist_fraction", "PET_total", "PET_intra", "PET_inter", "PET_intragt1000", "PET_intrale1000", "RNApANN", "RNADF.classifications", "singleS.Score", "singleG2M.Score", "singlePhase", "Age", "Age2", "RNA_snn_res.0.3", "seurat_clusters", "ATAC_snn_res.0.5", "nCount_PET", "nFeature_PET", "PET_snn_res.0.5", "RNA.weight", "ATAC.weight", "wsnn_res.0.3", "PET.weight", "wsnn3_res.0.3", "celltype2", "bigClass1", "bigClass2", "bigClass3", "totalATAC", "totalRNA", "umipercent", "atacpercent", "atacpercentcg", "intrapercent", "interpercent")) {P95[[i]] = NULL}
P95$groupid = "P95"
P95 <- NormalizeData(object = P95)
P95 <- FindVariableFeatures(object = P95, nfeatures = 3000)
P95 <- ScaleData(P95, features=rownames(P95))
P95 <- RunPCA(P95, features = VariableFeatures(object = P95))

for (i in c("X", "size", "cluster_color", "doublet.score", "cluster_id", "total.reads", "nonconf_mapped_reads", "unmapped_reads", "mapped_reads", "exp_component_name", "method", "Live_Cells", "Donor", "Total_Cells", "Live_percent", "Saturation", "Median_UMI_perCell", "Median_Genes_perCell", "gene.counts", "umi.counts", "Mean_Reads_perCell", "Lib_Cells", "Cell_Capture", "Lib_PassFail", "Lib_PCR_cycles", "Replicate_Lib", "Amp_PCR_cyles", "Amp_Date", "Amp_Name", "Gender", "Lib_type", "Region", "Seq_batch", "tube_barcode", "library_id", "aggr_num", "Lib_Date", "Lib_Name")) {SCG[[i]] = NULL}
SCG$groupid = "allen"


cca.peak.anchor <- FindIntegrationAnchors(object.list = c(P95, SCG), anchor.features = 3000, reduction = "cca")
saveRDS(cca.peak.anchor, "cca.peak.anchor.rds")
cca.Integrate <- IntegrateData(anchorset = cca.peak.anchor)
cca.Integrate <- ScaleData(cca.Integrate, features=rownames(cca.Integrate))
cca.Integrate <- RunPCA(cca.Integrate, features = VariableFeatures(cca.Integrate))
cca.Integrate <- RunTSNE(cca.Integrate, dims = 1:30, reduction="pca")
cca.Integrate <- RunUMAP(cca.Integrate, dims = 1:30, reduction="pca")

cca.Integrate$cluster_label[is.na(cca.Integrate$cluster_label)] = cca.Integrate$celltype[!is.na(cca.Integrate$celltype)]
cca.Integrate$subclass_label[is.na(cca.Integrate$subclass_label)] = cca.Integrate$celltype[!is.na(cca.Integrate$celltype)]
cca.Integrate$class_label[is.na(cca.Integrate$class_label)] = cca.Integrate$celltype[!is.na(cca.Integrate$celltype)]
pdf("plot.cca.peak.pdf", height=10, width=20)
p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group.by = "groupid") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group.by = "groupid") + ggtitle("umap")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, split.by = "groupid") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, split.by = "groupid") + ggtitle("umap")
print(p1)


p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group.by = "cluster_label", split.by = "groupid") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group.by = "cluster_label", split.by = "groupid") + ggtitle("umap")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group.by = "subclass_label", split.by = "groupid") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group.by = "subclass_label", split.by = "groupid") + ggtitle("umap")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group.by = "class_label", split.by = "groupid") + ggtitle("tsne")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group.by = "class_label", split.by = "groupid") + ggtitle("umap")
print(p1)
dev.off()
saveRDS(cca.Integrate, "cca.peak.anchor.integrate.rds")


cca.transfer.anchors <- FindTransferAnchors(reference = SCG, query = P95, dims = 1:30, reduction = 'cca' )
saveRDS(cca.transfer.anchors, "cca.transfer.anchors.rds")
cca.predicted.labels <- TransferData(anchorset = cca.transfer.anchors, refdata = SCG$subclass_label, weight.reduction = P95[['pca']], dims = 1:30)
P95 <- AddMetaData(object = P95, metadata = cca.predicted.labels)
P95 <- AddMetaData(object = P95, metadata = P95$predicted.id, col.name = 'subclass_label')
write.table(table(P95@meta.data[,c("subclass_label","celltype")]), file="cca.transfer.txt",sep="\t")
saveRDS(P95, file="cca.transfer.all.rds")
