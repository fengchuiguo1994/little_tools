library(Seurat)
library(Signac)
library(harmony)
library(ggplot2)
set.seed(1234)

RenameGenesSeurat_v2 <- function(obj,newnames,gene.use=NULL,de.assay="RNA") {
    print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
    lassays <- Assays(obj)
    assay.use <- obj@reductions$pca@assay.used
    DefaultAssay(obj) <- de.assay
    if (is.null(gene.use)) {
        all_genenames <- rownames(obj)
    }else{
        all_genenames <- gene.use
        obj <- subset(obj,features=gene.use)
    }

    order_name <- function(v1,v2,ref){
        v2 <- make.names(v2,unique=T)
        df1 <- data.frame(v1,v2)
        rownames(df1) <- df1$v1
        df1 <- df1[ref,]
        return(df1)
    }
    
    df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
    all_genenames <- df1$v1
    newnames <- df1$v2

    if ('SCT' %in% lassays) {
        if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
            obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
            rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
        }
    }
    change_assay <- function(a1=de.assay,obj,newnames=NULL,all_genenames=NULL){
        RNA <- obj@assays[a1][[1]]
        if (nrow(RNA) == length(newnames)) {
            if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
            if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
            if (length(RNA@var.features)) {
                df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
                all_genenames1 <- df1$v1
                newnames1 <- df1$v2
                RNA@var.features            <- newnames1
            }
            
            if (length(RNA@scale.data)){
                df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
                all_genenames1 <- df1$v1
                newnames1 <- df1$v2
                rownames(RNA@scale.data)    <- newnames1
            }
        } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
        
        obj@assays[a1][[1]] <- RNA
        return(obj)
    }

    for (a in lassays) {
        DefaultAssay(obj) <- a
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
    }

    hvg <- VariableFeatures(obj,assay=assay.use)
    if (length(obj@reductions$pca)){
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        rownames(obj@reductions$pca@feature.loadings) <- newnames1
    }
    return(obj)
}

RenameGenesSeurat <- function(obj,newnames,gene.use=NULL,de.assay="Spatial") {
    # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration.
    # It only changes obj@assays$RNA@counts, @data and @scale.data.
    print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
    lassays <- Assays(obj)
    #names(obj@assays)
    assay.use <- obj@reductions$pca@assay.used
    DefaultAssay(obj) <- de.assay
    if (is.null(gene.use)) {
        all_genenames <- rownames(obj)
    }else{
        all_genenames <- gene.use
        obj <- subset(obj,features=gene.use)
    }

    order_name <- function(v1,v2,ref){
        v2 <- make.names(v2,unique=T)
        df1 <- data.frame(v1,v2)
        rownames(df1) <- df1$v1
        df1 <- df1[ref,]
        return(df1)
    }

    df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
    all_genenames <- df1$v1
    newnames <- df1$v2

    if ('SCT' %in% lassays) {
        if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
            obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
            rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
        }
    }

    change_assay <- function(a1=de.assay,obj,newnames=NULL,all_genenames=NULL){
        RNA <- obj@assays[a1][[1]]
        if (nrow(RNA) == length(newnames)) {
            if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
            if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
            if (length(RNA@var.features)) {
                df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
                all_genenames1 <- df1$v1
                newnames1 <- df1$v2
                RNA@var.features            <- newnames1
            }
            if (length(RNA@scale.data)){
                df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
                all_genenames1 <- df1$v1
                newnames1 <- df1$v2
                rownames(RNA@scale.data)    <- newnames1
            }
        } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
        obj@assays[a1][[1]] <- RNA
        return(obj)
    }

    for (a in lassays) {
        DefaultAssay(obj) <- a
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
    }
    hvg <- VariableFeatures(obj,assay=assay.use)
    if (length(obj@reductions$pca)){
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
        df1 <- df1[rownames(obj@reductions$pca@feature.loadings),]
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        rownames(obj@reductions$pca@feature.loadings) <- newnames1
    }
    try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]])))
    return(obj)
}

SCG = readRDS("/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/8.mergeSampleBrainFinal/mouseBrain.P95.rds")
aa = SCG@meta.data
non = aa[aa$bigClass1!="Excitatory_Neurons" & aa$bigClass1!="Inhibitory_Neurons",]
neuron = aa[aa$bigClass1=="Excitatory_Neurons" | aa$bigClass1=="Inhibitory_Neurons",]
corneu = neuron[neuron$celltype2 == "TEGLU" | neuron$celltype2 == "TEINH" | neuron$celltype2 == "MSN",]
SCG$cellid = rownames(aa)
SCGflt = subset(SCG, subset=(cellid %in% c(rownames(non), rownames(corneu))))
SCG = SCGflt
mydat = as.data.frame(table(SCG@meta.data$celltype))
flt = as.character(mydat[mydat$Freq >= 50,]$Var1)
SCG = subset(SCG, subset = (celltype %in% flt))
SCG = subset(SCG, subset = (celltype != "VSMCA"))

SCG@reductions$pca = NULL
SCG@reductions$RNAPCA = NULL
SCG@reductions$tsne.rna = NULL
SCG@reductions$umap.rna = NULL
SCG@reductions$lsi = NULL
SCG@reductions$umap.atac = NULL
SCG@reductions$umap.pet = NULL
SCG@reductions$tsne.atac = NULL
SCG@reductions$tsne.pet = NULL
SCG@reductions$wnn.umap = NULL
SCG@reductions$wnn.tsne = NULL
SCG@reductions$wnn3.umap = NULL
SCG@reductions$wnn3.tsne = NULL
SCG@reductions$petlsi = NULL


mydat = readRDS("c92c5b85-61d4-4228-b8fc-24ce4a757b59.rds")
mydat

pdf("mouseBrain.metadata.raw.plot.pdf", height=10, width=20)
p1 = DimPlot(mydat, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("pca")
print(p1)
p1 = DimPlot(mydat, reduction ="pca_harmony", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("pca_harmony")
print(p1)
p1 = DimPlot(mydat, reduction ="pca_orig", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("pca_orig")
print(p1)
p1 = DimPlot(mydat, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("umap")
print(p1)
p1 = DimPlot(mydat, reduction ="spatial_coords", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("spatial_coords")
print(p1)
dev.off()
write.table(mydat@meta.data, file="mouseBrain.metadata.txt", quote=F, sep="\t")

oldname = rownames(mydat@assays$RNA@data)
mydat@assays$RNA@data = as(mydat@assays$RNA@data,"dgCMatrix")
maplist = read.table("/data/home/ruanlab/huangxingyu/Tools/cellranger-7.1.0/genome/refdata-gex-mm10-2020-A/genes/mm10.info",header=F)
names(maplist) = c("chromosome", "start", "end", "geneid", "score", "strand", "genename", "type")
rownames(maplist) = maplist$geneid
maplistflt = maplist[oldname,]

mydat2 <- RenameGenesSeurat_v2(obj=mydat,newnames=maplist$genename,gene.use=maplist$geneid)
sharegene = intersect(rownames(mydat2@assays$RNA@counts),rownames(SCG@assays$RNA@counts))

SCGflt = SCG[sharegene,]
mydatflt = mydat2[sharegene,]

SCGflt <- NormalizeData(object = SCGflt)
mydatflt <- NormalizeData(object = mydatflt)
SCGflt <- FindVariableFeatures(object = SCGflt, nfeatures = 3000)
mydatflt[["RNA"]]@meta.features <- data.frame(row.names = rownames(mydatflt[["RNA"]])) # for solve error # https://github.com/satijalab/seurat/issues/2317
mydatflt <- FindVariableFeatures(object = mydatflt, nfeatures = 3000)
SCGflt <- ScaleData(SCGflt, features=rownames(SCGflt))
mydatflt <- ScaleData(mydatflt, features=rownames(mydatflt))

mydatflt@reductions$pca = NULL
mydatflt@reductions$pca_harmony = NULL
mydatflt@reductions$pca_orig = NULL
mydatflt@reductions$umap = NULL
mydatflt@reductions$spatial_coords = NULL

SCGflt <- RunPCA(SCGflt, features = VariableFeatures(object = SCGflt))
mydatflt <- RunPCA(mydatflt, features = VariableFeatures(object = mydatflt))
SCGflt <- RunPCA(SCGflt, features = VariableFeatures(object = SCGflt), reduction.name="pcaraw", reduction.key = "pcarawPC_")
mydatflt <- RunPCA(mydatflt, features = VariableFeatures(object = mydatflt), reduction.name="pcaraw", reduction.key = "pcarawPC_")
SCGflt <- RunTSNE(SCGflt, dims = 1:30)
mydatflt <- RunTSNE(mydatflt, dims = 1:30)
SCGflt <- RunUMAP(SCGflt, dims = 1:30)
mydatflt <- RunUMAP(mydatflt, dims = 1:30)
SCGflt <- RunTSNE(SCGflt, dims = 1:30, reduction="pca", reduction.name="tsneraw", reduction.key = "tsnerawTSNE_")
mydatflt <- RunTSNE(mydatflt, dims = 1:30, reduction="pca", reduction.name="tsneraw", reduction.key = "tsnerawTSNE_")
SCGflt <- RunUMAP(SCGflt, dims = 1:30, reduction="pca", reduction.name="umapraw", reduction.key = "tsnerawUMAP_")
mydatflt <- RunUMAP(mydatflt, dims = 1:30, reduction="pca", reduction.name="umapraw", reduction.key = "tsnerawUMAP_")

SCGflt <- FindNeighbors(SCGflt, dims = 1:30)
mydatflt <- FindNeighbors(mydatflt, dims = 1:30)
SCGflt <- FindClusters(SCGflt, resolution = 0.8)
mydatflt <- FindClusters(mydatflt, resolution = 0.8)

pdf("mouseBrain.zhuang.plot.pdf", height=10, width=20)
p1 = DimPlot(mydatflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("pca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("tsne")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("umap")
print(p1)

p1 = DimPlot(mydatflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "donor_id") + ggtitle("pca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "donor_id") + ggtitle("tsne")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "donor_id") + ggtitle("umap")
print(p1)

p1 = DimPlot(mydatflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("pca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("tsne")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("umap")
print(p1)
dev.off()

pdf("mouseBrain.SCG.plot.pdf", height=10, width=20)
p1 = DimPlot(SCGflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "bigClass1") + ggtitle("pca")
print(p1)
p1 = DimPlot(SCGflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "bigClass1") + ggtitle("tsne")
print(p1)
p1 = DimPlot(SCGflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "bigClass1") + ggtitle("umap")
print(p1)

p1 = DimPlot(SCGflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "orig.ident") + ggtitle("pca")
print(p1)
p1 = DimPlot(SCGflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "orig.ident") + ggtitle("tsne")
print(p1)
p1 = DimPlot(SCGflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "orig.ident") + ggtitle("umap")
print(p1)

p1 = DimPlot(SCGflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("pca")
print(p1)
p1 = DimPlot(SCGflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("tsne")
print(p1)
p1 = DimPlot(SCGflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("umap")
print(p1)
dev.off()


pdf("mouseBrain.zhuang.plot.split.pdf", height=10, width=40)
p1 = DimPlot(mydatflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot", split.by = "donor_id") + ggtitle("pca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot", split.by = "donor_id") + ggtitle("tsne")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot", split.by = "donor_id") + ggtitle("umap")
print(p1)

p1 = DimPlot(mydatflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "donor_id", split.by = "donor_id") + ggtitle("pca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "donor_id", split.by = "donor_id") + ggtitle("tsne")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "donor_id", split.by = "donor_id") + ggtitle("umap")
print(p1)

p1 = DimPlot(mydatflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "donor_id") + ggtitle("pca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "donor_id") + ggtitle("tsne")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "donor_id") + ggtitle("umap")
print(p1)
dev.off()

pdf("mouseBrain.SCG.plot.split.pdf", height=10, width=40)
p1 = DimPlot(SCGflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "bigClass1", split.by = "orig.ident") + ggtitle("pca")
print(p1)
p1 = DimPlot(SCGflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "bigClass1", split.by = "orig.ident") + ggtitle("tsne")
print(p1)
p1 = DimPlot(SCGflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "bigClass1", split.by = "orig.ident") + ggtitle("umap")
print(p1)

p1 = DimPlot(SCGflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "orig.ident", split.by = "orig.ident") + ggtitle("pca")
print(p1)
p1 = DimPlot(SCGflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "orig.ident", split.by = "orig.ident") + ggtitle("tsne")
print(p1)
p1 = DimPlot(SCGflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "orig.ident", split.by = "orig.ident") + ggtitle("umap")
print(p1)

p1 = DimPlot(SCGflt, reduction ="pca", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "orig.ident") + ggtitle("pca")
print(p1)
p1 = DimPlot(SCGflt, reduction ="tsne", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "orig.ident") + ggtitle("tsne")
print(p1)
p1 = DimPlot(SCGflt, reduction ="umap", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "orig.ident") + ggtitle("umap")
print(p1)
dev.off()


mydatflt <- RunHarmony(mydatflt, group.by.vars = "age", reduction = "pcaraw", reduction.save = "rnaharmonypca", max.iter.harmony = 100)
# mydatflt <- RunHarmony(mydatflt, group.by.vars = "age", reduction = "pcaraw", reduction.save = "pca", max.iter.harmony = 100)
mydatflt@reductions$pca = mydatflt@reductions$rnaharmonypca
mydatflt <- RunTSNE(mydatflt, dims = 1:30, reduction="pca", reduction.name="tsneharmony", reduction.key = "tsneharmonyTSNE_")
mydatflt <- RunUMAP(mydatflt, dims = 1:30, reduction="pca", reduction.name="umapharmony", reduction.key = "tsneharmonyUMAP_")

pdf("mouseBrain.zhuang.plot.harmony.pdf", height=10, width=20)
p1 = DimPlot(mydatflt, reduction ="rnaharmonypca", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsneharmony", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("tsneharmony")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umapharmony", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot") + ggtitle("umapharmony")
print(p1)

p1 = DimPlot(mydatflt, reduction ="rnaharmonypca", label = TRUE, label.size=6, raster=FALSE, group = "donor_id") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsneharmony", label = TRUE, label.size=6, raster=FALSE, group = "donor_id") + ggtitle("tsneharmony")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umapharmony", label = TRUE, label.size=6, raster=FALSE, group = "donor_id") + ggtitle("umapharmony")
print(p1)

p1 = DimPlot(mydatflt, reduction ="rnaharmonypca", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsneharmony", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("tsneharmony")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umapharmony", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters") + ggtitle("umapharmony")
print(p1)
dev.off()

pdf("mouseBrain.zhuang.plot.harmony.split.pdf", height=10, width=40)
p1 = DimPlot(mydatflt, reduction ="rnaharmonypca", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot", split.by = "donor_id") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsneharmony", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot", split.by = "donor_id") + ggtitle("tsneharmony")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umapharmony", label = TRUE, label.size=6, raster=FALSE, group = "clust_annot", split.by = "donor_id") + ggtitle("umapharmony")
print(p1)

p1 = DimPlot(mydatflt, reduction ="rnaharmonypca", label = TRUE, label.size=6, raster=FALSE, group = "donor_id", split.by = "donor_id") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsneharmony", label = TRUE, label.size=6, raster=FALSE, group = "donor_id", split.by = "donor_id") + ggtitle("tsneharmony")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umapharmony", label = TRUE, label.size=6, raster=FALSE, group = "donor_id", split.by = "donor_id") + ggtitle("umapharmony")
print(p1)

p1 = DimPlot(mydatflt, reduction ="rnaharmonypca", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "donor_id") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(mydatflt, reduction ="tsneharmony", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "donor_id") + ggtitle("tsneharmony")
print(p1)
p1 = DimPlot(mydatflt, reduction ="umapharmony", label = TRUE, label.size=6, raster=FALSE, group = "seurat_clusters", split.by = "donor_id") + ggtitle("umapharmony")
print(p1)
dev.off()

saveRDS(mydatflt, file="zhuang.cg.rds")
saveRDS(SCGflt, file="P95.cg.rds")
# mydatflt = readRDS(file="zhuang.cg.rds")
# SCGflt = readRDS(file="P95.cg.rds")

cca.peak.anchor <- FindIntegrationAnchors(object.list = c(mydatflt, SCGflt), anchor.features = 3000, reduction = "cca")
saveRDS(cca.peak.anchor, "cca.peak.anchor.rds")
cca.Integrate <- IntegrateData(anchorset = cca.peak.anchor)
cca.Integrate <- ScaleData(cca.Integrate, features=rownames(cca.Integrate))
cca.Integrate <- RunPCA(cca.Integrate, features = VariableFeatures(cca.Integrate))
cca.Integrate <- RunPCA(cca.Integrate, features = VariableFeatures(cca.Integrate), reduction.name="pcaraw", reduction.key = "pcarawPC_")
cca.Integrate <- RunTSNE(cca.Integrate, dims = 1:30, reduction="pca")
cca.Integrate <- RunUMAP(cca.Integrate, dims = 1:30, reduction="pca")
cca.Integrate <- RunTSNE(cca.Integrate, dims = 1:30, reduction="pca", reduction.name="tsneraw", reduction.key = "tsnerawTSNE_")
cca.Integrate <- RunUMAP(cca.Integrate, dims = 1:30, reduction="pca", reduction.name="umapraw", reduction.key = "tsnerawUMAP_")
aaa = cca.Integrate@meta.data
aaa$combineanno = "#"
aaa[!is.na(aaa$celltype),]$combineanno = aaa[!is.na(aaa$celltype),]$celltype
aaa[!is.na(aaa$clust_annot),]$combineanno = aaa[!is.na(aaa$clust_annot),]$clust_annot
aaa$combinegroup = "#"
aaa[!is.na(aaa$age),]$combinegroup = "zhuang"
aaa[is.na(aaa$age),]$combinegroup = "P95"
cca.Integrate$combineanno = aaa$combineanno
cca.Integrate$combinegroup = aaa$combinegroup
pdf("mouseBrain.zhuang.plot.cca.peak.pdf", height=10, width=20)
p1 = DimPlot(cca.Integrate, reduction ="pcaraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsneraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("tsneraw")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umapraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("umapraw")
print(p1)
dev.off()
pdf("mouseBrain.zhuang.plot.cca.peak.split.pdf", height=10, width=40)
p1 = DimPlot(cca.Integrate, reduction ="pcaraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="tsneraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("tsneraw")
print(p1)
p1 = DimPlot(cca.Integrate, reduction ="umapraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("umapraw")
print(p1)
dev.off()
saveRDS(cca.Integrate, file = "cca.Integrate.all.rds")


rpca.peak.anchor <- FindIntegrationAnchors(object.list = c(mydatflt, SCGflt), anchor.features = 3000, reduction = "rpca")
saveRDS(rpca.peak.anchor, "rpca.peak.anchor.rds")
rpca.Integrate <- IntegrateData(anchorset = rpca.peak.anchor)
rpca.Integrate <- ScaleData(rpca.Integrate, features=rownames(rpca.Integrate))
rpca.Integrate <- RunPCA(rpca.Integrate, features = VariableFeatures(rpca.Integrate))
rpca.Integrate <- RunPCA(rpca.Integrate, features = VariableFeatures(rpca.Integrate), reduction.name="pcaraw", reduction.key = "pcarawPC_")
rpca.Integrate <- RunTSNE(rpca.Integrate, dims = 1:30, reduction="pca")
rpca.Integrate <- RunUMAP(rpca.Integrate, dims = 1:30, reduction="pca")
rpca.Integrate <- RunTSNE(rpca.Integrate, dims = 1:30, reduction="pca", reduction.name="tsneraw", reduction.key = "tsnerawTSNE_")
rpca.Integrate <- RunUMAP(rpca.Integrate, dims = 1:30, reduction="pca", reduction.name="umapraw", reduction.key = "tsnerawUMAP_")
aaa = rpca.Integrate@meta.data
aaa$combineanno = "#"
aaa[!is.na(aaa$celltype),]$combineanno = aaa[!is.na(aaa$celltype),]$celltype
aaa[!is.na(aaa$clust_annot),]$combineanno = aaa[!is.na(aaa$clust_annot),]$clust_annot
aaa$combinegroup = "#"
aaa[!is.na(aaa$age),]$combinegroup = "zhuang"
aaa[is.na(aaa$age),]$combinegroup = "P95"
rpca.Integrate$combineanno = aaa$combineanno
rpca.Integrate$combinegroup = aaa$combinegroup
pdf("mouseBrain.zhuang.plot.rpca.peak.pdf", height=10, width=20)
p1 = DimPlot(rpca.Integrate, reduction ="pcaraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(rpca.Integrate, reduction ="tsneraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("tsneraw")
print(p1)
p1 = DimPlot(rpca.Integrate, reduction ="umapraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("umapraw")
print(p1)
dev.off()
pdf("mouseBrain.zhuang.plot.rpca.peak.split.pdf", height=10, width=40)
p1 = DimPlot(rpca.Integrate, reduction ="pcaraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(rpca.Integrate, reduction ="tsneraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("tsneraw")
print(p1)
p1 = DimPlot(rpca.Integrate, reduction ="umapraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("umapraw")
print(p1)
dev.off()
saveRDS(rpca.Integrate, file = "rpca.Integrate.all.rds")


harmonymerge <- merge(mydatflt, y = c(SCGflt), project = "mouseBrain")
harmonymerge <- NormalizeData(object = harmonymerge)
harmonymerge <- FindVariableFeatures(object = harmonymerge, nfeatures = 3000)
harmonymerge <- ScaleData(harmonymerge, features=rownames(harmonymerge))
harmonymerge <- RunPCA(harmonymerge, features = VariableFeatures(object = harmonymerge))
harmonymerge <- RunPCA(harmonymerge, features = VariableFeatures(object = harmonymerge), reduction.name="pcaraw", reduction.key = "pcarawPC_")
harmonymerge <- RunTSNE(harmonymerge, dims = 1:30)
harmonymerge <- RunUMAP(harmonymerge, dims = 1:30)
harmonymerge <- RunTSNE(harmonymerge, dims = 1:30, reduction="pca", reduction.name="tsneraw", reduction.key = "tsnerawTSNE_")
harmonymerge <- RunUMAP(harmonymerge, dims = 1:30, reduction="pca", reduction.name="umapraw", reduction.key = "tsnerawUMAP_")

aaa = harmonymerge@meta.data
aaa$combineanno = "#"
aaa[!is.na(aaa$celltype),]$combineanno = aaa[!is.na(aaa$celltype),]$celltype
aaa[!is.na(aaa$clust_annot),]$combineanno = aaa[!is.na(aaa$clust_annot),]$clust_annot
aaa$combinegroup = "#"
aaa[!is.na(aaa$age),]$combinegroup = "zhuang"
aaa[is.na(aaa$age),]$combinegroup = "P95"
aaa$agegroup = "#"
aaa[!is.na(aaa$age),]$agegroup = aaa[!is.na(aaa$age),]$age
aaa[is.na(aaa$age),]$agegroup = aaa[is.na(aaa$age),]$Age
harmonymerge$combineanno = aaa$combineanno
harmonymerge$combinegroup = aaa$combinegroup
harmonymerge$agegroup = aaa$agegroup

harmonymerge <- RunHarmony(harmonymerge, group.by.vars = "agegroup", reduction = "pcaraw", reduction.save = "rnaharmonypca", max.iter.harmony = 100)
harmonymerge@reductions$pca = harmonymerge@reductions$rnaharmonypca
harmonymerge <- RunTSNE(harmonymerge, dims = 1:30, reduction="rnaharmonypca", reduction.name="tsneharmony", reduction.key = "tsneharmonyTSNE_")
harmonymerge <- RunUMAP(harmonymerge, dims = 1:30, reduction="rnaharmonypca", reduction.name="umapharmony", reduction.key = "tsneharmonyUMAP_")

pdf("mouseBrain.zhuang.plot.harmonymerge.pdf", height=10, width=20)
p1 = DimPlot(harmonymerge, reduction ="pcaraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("pcaraw")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="tsneraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("tsneraw")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="umapraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("umapraw")
print(p1)

p1 = DimPlot(harmonymerge, reduction ="rnaharmonypca", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="tsneharmony", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("tsneharmony")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="umapharmony", label = TRUE, label.size=6, raster=FALSE, group = "combineanno") + ggtitle("umapharmony")
print(p1)
dev.off()
pdf("mouseBrain.zhuang.plot.harmonymerge.split.pdf", height=10, width=40)
p1 = DimPlot(harmonymerge, reduction ="pcaraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("pcaraw")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="tsneraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("tsneraw")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="umapraw", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("umapraw")
print(p1)

p1 = DimPlot(harmonymerge, reduction ="rnaharmonypca", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("rnaharmonypca")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="tsneharmony", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("tsneharmony")
print(p1)
p1 = DimPlot(harmonymerge, reduction ="umapharmony", label = TRUE, label.size=6, raster=FALSE, group = "combineanno", split.by = "combinegroup") + ggtitle("umapharmony")
print(p1)
dev.off()
saveRDS(harmonymerge, file = "harmonymerge.all.rds")


pcaproject.transfer.anchors <- FindTransferAnchors(reference = mydatflt, query = SCGflt, dims = 1:30, reduction = 'pcaproject' )
saveRDS(pcaproject.transfer.anchors, "pcaproject.transfer.anchors.rds")
pcaproject.predicted.labels <- TransferData(anchorset = pcaproject.transfer.anchors, refdata = mydatflt$clust_annot, weight.reduction = SCGflt[['pca']], dims = 1:30)
SCGflt <- AddMetaData(object = SCGflt, metadata = pcaproject.predicted.labels)
SCGflt <- AddMetaData(object = SCGflt, metadata = SCGflt$predicted.id, col.name = 'clust_annot')
saveRDS(SCGflt, file="pcaproject.transfer.rds")
write.table(table(SCGflt@meta.data[,c("clust_annot","celltype")]), file="pcaproject.transfer.txt",sep="\t")
saveRDS(SCGflt, file="pcaproject.transfer.all.rds")


rpca.transfer.anchors <- FindTransferAnchors(reference = mydatflt, query = SCGflt, dims = 1:30, reduction = 'rpca' )
saveRDS(rpca.transfer.anchors, "rpca.transfer.anchors.rds")
rpca.predicted.labels <- TransferData(anchorset = rpca.transfer.anchors, refdata = mydatflt$clust_annot, weight.reduction = SCGflt[['pca']], dims = 1:30)
SCGflt <- AddMetaData(object = SCGflt, metadata = rpca.predicted.labels)
SCGflt <- AddMetaData(object = SCGflt, metadata = SCGflt$predicted.id, col.name = 'clust_annot')
saveRDS(SCGflt, file="rpca.transfer.rds")
write.table(table(SCGflt@meta.data[,c("clust_annot","celltype")]), file="rpca.transfer.txt",sep="\t")
saveRDS(SCGflt, file="rpca.transfer.all.rds")


cca.transfer.anchors <- FindTransferAnchors(reference = mydatflt, query = SCGflt, dims = 1:30, reduction = 'cca' )
saveRDS(cca.transfer.anchors, "cca.transfer.anchors.rds")
cca.predicted.labels <- TransferData(anchorset = cca.transfer.anchors, refdata = mydatflt$clust_annot, weight.reduction = SCGflt[['pca']], dims = 1:30)
SCGflt <- AddMetaData(object = SCGflt, metadata = cca.predicted.labels)
SCGflt <- AddMetaData(object = SCGflt, metadata = SCGflt$predicted.id, col.name = 'clust_annot')
saveRDS(SCGflt, file="cca.transfer.rds")
write.table(table(SCGflt@meta.data[,c("clust_annot","celltype")]), file="cca.transfer.txt",sep="\t")
saveRDS(SCGflt, file="cca.transfer.all.rds")