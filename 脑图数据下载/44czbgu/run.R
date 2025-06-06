library(Seurat)
library(dplyr)
library(ggplot2)
library(Signac)
set.seed(1234)

counts <- Read10X(data.dir = "44czbgu")
scRNA <- CreateSeuratObject(counts = counts, project = "44czbgu", assay = "RNA", min.cells = 3, min.features = 3)
saveRDS(scRNA, file = "44czbgu.raw.data.rds")

metadata = read.csv("sample_metadata.csv",header=T)
rownames(metadata) = metadata$sample_name
metadata = metadata[rownames(scRNA@meta.data),]
# table(rownames(scRNA@meta.data) == rownames(metadata))
metadata = metadata[,-1]
scRNA = AddMetaData(scRNA, metadata)

anno1 = read.csv("cluster.membership.csv", header=T)
names(anno1)[2] = "cluster_id"
anno2 = read.csv("cluster.annotation.csv", header=T)
anno = left_join(anno1,anno2,by="cluster_id")

fltlist = intersect(metadata$sample_name, anno$X)
droplist = setdiff(metadata$sample_name, anno$X)

scRNAflt = subset(scRNA, subset=(sample_name %in% fltlist))
scRNAdrop = subset(scRNA, subset=(sample_name %in% droplist))
rownames(anno) = anno$X
anno = anno[rownames(scRNAflt@meta.data),]
# table(rownames(scRNAflt@meta.data) == rownames(anno))
annod = anno[,-1]
scRNAflt = AddMetaData(scRNAflt, annod)


scRNAflt <- NormalizeData(object = scRNAflt)
scRNAflt <- FindVariableFeatures(object = scRNAflt, nfeatures = 3000)
scRNAflt <- ScaleData(scRNAflt, features=rownames(scRNAflt))
scRNAflt <- RunPCA(scRNAflt, features = VariableFeatures(object = scRNAflt))
scRNAflt <- RunUMAP(scRNAflt, dims = 1:30)
scRNAflt <- RunTSNE(scRNAflt, dims = 1:30)

pdf("44czbgu.class.pdf", height=10, width=20)
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

saveRDS(scRNAflt, file = "44czbgu.cluster.rds")




scRNAflt = readRDS(file = "44czbgu.cluster.rds")
# aaa = scRNAflt@meta.data
# head(aaa[aaa$subclass_label=="doublet" | aaa$subclass_label=="Low Quality" | aaa$subclass_label=="Batch" | aaa$class_label=="Low Quality",])
# table(aaa[aaa$subclass_label=="doublet",]$cluster_label)
# table(aaa[aaa$subclass_label=="Batch",]$cluster_label)
# table(aaa[aaa$subclass_label=="Low Quality",]$cluster_label)
# table(aaa[aaa$class_label=="Low Quality",]$cluster_label)
SCG = subset(scRNAflt, subset=(subclass_label!="doublet" & subclass_label!="Low Quality" & subclass_label!="Batch" & class_label!="Low Quality"))
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

pdf("44czbgu.final.class.pdf", height=10, width=20)
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
saveRDS(SCG, file = "44czbgu.final.cluster.rds")

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

for (i in c("sample_name", "nUMI", "nGene", "dataset", "QC", "cluster", "Allen.cluster_id", "Allen.cluster_color", "comb.QC", "cluster_id", "cluster_color", "size", "gene.counts", "umi.counts", "Broad.QC.doublet", "Broad.QC.Mito", "Broad.passQC", "MALE", "Comb.QC", "cl", "Allen.cluster_label", "Allen.class_label", "Allen.subclass_label")) {SCG[[i]] = NULL}
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
write.table(P95@meta.data, file="cca.transfer.metadata.txt",sep="\t")
saveRDS(P95, file="cca.transfer.all.rds")










#####################################################
## reintegrate again
SCG = readRDS(file = "44czbgu.final.cluster.rds")
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
P95 = subset(P95, subset = (celltype %in% c("ACTE1","ACTE2","OPC","MFOL2","MOL1","MOL2","MOL3","MGL1","MGL2","MGL3","PER1","PER3","PVM1","TEGLU1","TEGLU2","TEGLU3","TEGLU4","TEGLU5","TEGLU6","TEGLU7","TEGLU8","TEGLU10","TEGLU11","TEGLU12","TEINH4","TEINH15","TEINH17","TEINH18","TEINH19","VECV","VLMC1","VLMC2")))

P95 = DietSeurat(P95, counts = TRUE, scale.data = F, assays = "RNA")
for (i in c("percent.mt", "percent.ribo", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "TSS.percentile", "nucleosome_signal", "nucleosome_percentile", "blacklist_fraction", "PET_total", "PET_intra", "PET_inter", "PET_intragt1000", "PET_intrale1000", "RNApANN", "RNADF.classifications", "singleS.Score", "singleG2M.Score", "singlePhase", "Age", "Age2", "RNA_snn_res.0.3", "seurat_clusters", "ATAC_snn_res.0.5", "nCount_PET", "nFeature_PET", "PET_snn_res.0.5", "RNA.weight", "ATAC.weight", "wsnn_res.0.3", "PET.weight", "wsnn3_res.0.3", "celltype2", "bigClass1", "bigClass2", "bigClass3", "totalATAC", "totalRNA", "umipercent", "atacpercent", "atacpercentcg", "intrapercent", "interpercent")) {P95[[i]] = NULL}
P95$groupid = "P95"
P95 <- NormalizeData(object = P95)
P95 <- FindVariableFeatures(object = P95, nfeatures = 3000)
P95 <- ScaleData(P95, features=rownames(P95))
P95 <- RunPCA(P95, features = VariableFeatures(object = P95))

for (i in c("sample_name", "nUMI", "nGene", "dataset", "QC", "cluster", "Allen.cluster_id", "Allen.cluster_color", "comb.QC", "cluster_id", "cluster_color", "size", "gene.counts", "umi.counts", "Broad.QC.doublet", "Broad.QC.Mito", "Broad.passQC", "MALE", "Comb.QC", "cl", "Allen.cluster_label", "Allen.class_label", "Allen.subclass_label")) {SCG[[i]] = NULL}
SCG$groupid = "allen"

cca.peak.anchor <- FindIntegrationAnchors(object.list = c(P95, SCG), anchor.features = 3000, reduction = "cca")
saveRDS(cca.peak.anchor, "cca.peak.anchor.reintegrate.rds")
cca.Integrate <- IntegrateData(anchorset = cca.peak.anchor)
cca.Integrate <- ScaleData(cca.Integrate, features=rownames(cca.Integrate))
cca.Integrate <- RunPCA(cca.Integrate, features = VariableFeatures(cca.Integrate))
cca.Integrate <- RunTSNE(cca.Integrate, dims = 1:30, reduction="pca")
cca.Integrate <- RunUMAP(cca.Integrate, dims = 1:30, reduction="pca")

cca.Integrate$cluster_label[is.na(cca.Integrate$cluster_label)] = cca.Integrate$celltype[!is.na(cca.Integrate$celltype)]
cca.Integrate$subclass_label[is.na(cca.Integrate$subclass_label)] = cca.Integrate$celltype[!is.na(cca.Integrate$celltype)]
cca.Integrate$class_label[is.na(cca.Integrate$class_label)] = cca.Integrate$celltype[!is.na(cca.Integrate$celltype)]
pdf("plot.cca.peak.reintegrate.pdf", height=10, width=20)
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
saveRDS(cca.Integrate, "cca.peak.anchor.integrate.reintegrate.rds")


cca.transfer.anchors <- FindTransferAnchors(reference = SCG, query = P95, dims = 1:30, reduction = 'cca' )
saveRDS(cca.transfer.anchors, "cca.transfer.anchors.reintegrate.rds")
cca.predicted.labels <- TransferData(anchorset = cca.transfer.anchors, refdata = SCG$subclass_label, weight.reduction = P95[['pca']], dims = 1:30)
P95 <- AddMetaData(object = P95, metadata = cca.predicted.labels)
P95 <- AddMetaData(object = P95, metadata = P95$predicted.id, col.name = 'subclass_label')
write.table(table(P95@meta.data[,c("subclass_label","celltype")]), file="cca.transfer.reintegrate.txt",sep="\t")
write.table(P95@meta.data, file="cca.transfer.metadata.reintegrate.txt",sep="\t")
saveRDS(P95, file="cca.transfer.all.reintegrate.rds")



P95 <- RunTSNE(P95, dims = 1:30, reduction="pca")
P95 <- RunUMAP(P95, dims = 1:30, reduction="pca")
pdf("plot.cca.peak.reintegrate.2.pdf", height=10, width=20)
p1 = DimPlot(P95, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group.by = "subclass_label", split.by = "groupid") + ggtitle("umap")
print(p1)
dev.off()
