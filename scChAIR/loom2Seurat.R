library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(clustree)

set.seed(1234)
pdf("MouseBrain.reference.pdf",width=15)

loom_mouse_brain_cell <- Connect(filename = "l5_all.loom",mode = "r")
loom_mouse_brain_cell.seurat <- as.Seurat(loom_mouse_brain_cell,cells = "ClusterName", features = "Gene")
saveRDS(loom_mouse_brain_cell.seurat, file="mouseBrain.reference.rds")
loom_mouse_brain_cell.seurat <- subset(loom_mouse_brain_cell.seurat, subset = !(Tissue == "DRG" | Tissue == "ENS" | Tissue == "SC" | Tissue == "Sympath" | Region == "Spinal cord") )
loom_mouse_brain_cell.seurat <- AddMetaData(loom_mouse_brain_cell.seurat, sapply(strsplit(rownames(loom_mouse_brain_cell.seurat@meta.data),".",fixed = TRUE), function(x) x[1]), col.name = "celltype")

x = "reference"

loom_mouse_brain_cell.seurat <- NormalizeData(object = loom_mouse_brain_cell.seurat, verbose=F)
loom_mouse_brain_cell.seurat <- FindVariableFeatures(object = loom_mouse_brain_cell.seurat,nfeatures = 3000, verbose=F)
loom_mouse_brain_cell.seurat <- ScaleData(loom_mouse_brain_cell.seurat, features=rownames(loom_mouse_brain_cell.seurat))
loom_mouse_brain_cell.seurat <- RunPCA(object = loom_mouse_brain_cell.seurat,features = VariableFeatures(loom_mouse_brain_cell.seurat), verbose = F)

loom_mouse_brain_cell.seurat <- JackStraw(loom_mouse_brain_cell.seurat, num.replicate = 100,dims = 50)
loom_mouse_brain_cell.seurat <- ScoreJackStraw(loom_mouse_brain_cell.seurat, dims = 1:50) 
p1 = JackStrawPlot(loom_mouse_brain_cell.seurat, dims = 1:50) + ggtitle(x)
p2 = ElbowPlot(loom_mouse_brain_cell.seurat, ndims = 50) + ggtitle(x)
print(p1 + p2)

loom_mouse_brain_cell.seurat <- FindNeighbors(loom_mouse_brain_cell.seurat, dims = 1:50)
loom_mouse_brain_cell.seurat <- FindClusters(loom_mouse_brain_cell.seurat, resolution = c(seq(0,2,.1)))
p1 = clustree(loom_mouse_brain_cell.seurat@meta.data, prefix = "RNA_snn_res.", node_colour_aggr = "median") + coord_flip()
print(p1)
for (x in grep("RNA_snn_res.",colnames(loom_mouse_brain_cell.seurat@meta.data),value = T)){loom_mouse_brain_cell.seurat[[x]] <- NULL}
loom_mouse_brain_cell.seurat <- FindClusters(loom_mouse_brain_cell.seurat, resolution = 0.3)

loom_mouse_brain_cell.seurat <- RunTSNE(loom_mouse_brain_cell.seurat, dims = 1:50, verbose = F, reduction.name="tsne.rna", reduction.key = "RNATSNE_") 
loom_mouse_brain_cell.seurat <- RunUMAP(loom_mouse_brain_cell.seurat, dims = 1:50, verbose = F, reduction.name="umap.rna", reduction.key = "RNAUMAP_")
p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=6, raster=FALSE) + ggtitle(paste(x,"rawRNAtsne",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=6, raster=FALSE) + ggtitle(paste(x,"rwaRNAumap",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=6, raster=FALSE, group.by="TaxonomyRank2") + ggtitle(paste(x,"rawRNAtsneTaxonomyRank2",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=6, raster=FALSE, group.by="TaxonomyRank2") + ggtitle(paste(x,"rwaRNAumapTaxonomyRank2",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=6, raster=FALSE, group.by="Taxonomy_group") + ggtitle(paste(x,"rawRNAtsneTaxonomy_group",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=6, raster=FALSE, group.by="Taxonomy_group") + ggtitle(paste(x,"rwaRNAumapTaxonomy_group",sep="."))
print(p4)


















p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomyRank1") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsneTaxonomyRank1",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomyRank1") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapTaxonomyRank1",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomyRank2") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsneTaxonomyRank2",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomyRank2") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapTaxonomyRank2",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomyRank3") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsneTaxonomyRank3",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomyRank3") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapTaxonomyRank3",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomyRank4") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsneTaxonomyRank4",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomyRank4") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapTaxonomyRank4",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomySymbol") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsneTaxonomySymbol",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="TaxonomySymbol") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapTaxonomySymbol",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="Taxonomy_group") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsneTaxonomy_group",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="Taxonomy_group") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapTaxonomy_group",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="Class") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsneClass",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="Class") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapClass",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="Subclass") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsneSubclass",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="Subclass") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapSubclass",sep="."))
print(p4)

p3 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="tsne.rna", label = TRUE, label.size=1, raster=FALSE, group.by="celltype") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rawRNAtsnecelltype",sep="."))
print(p3)
p4 = DimPlot(loom_mouse_brain_cell.seurat, reduction ="umap.rna", label = TRUE, label.size=1, raster=FALSE, group.by="celltype") + theme(legend.title = element_text(size=1), legend.text = element_text(size=6)) + ggtitle(paste(x,"rwaRNAumapcelltype",sep="."))
print(p4)

dev.off()

saveRDS(loom_mouse_brain_cell.seurat, file="mouseBrain.brainPart.reference.rds")