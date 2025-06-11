library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

atac_dir = "/data/home/ruanlab/huangxingyu/scChAIR/ChIAV3/example/scATAC/ATAC10XlocalInOutPutATAC/ATAC"
f_binary_mat = Matrix::readMM(paste0(atac_dir, "/matrix.mtx.gz"))
regions.names = read.delim(paste0(atac_dir, "/features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
cells.names = read.delim(paste0(atac_dir, "/barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
colnames(f_binary_mat) = cells.names$V1
rownames(f_binary_mat) = regions.names$V1
# f_binary_mat@x[f_binary_mat@x > 0] <- 1
message('Matrix size:\n', 'rows ', f_binary_mat@Dim[1], '\ncolumns ', f_binary_mat@Dim[2])

pdf("test.pdf")
n_cells_with_site = rowSums(f_binary_mat)
options(repr.plot.width = 8, repr.plot.height = 4)
par(mfrow = c(1, 2))
hist(log10(n_cells_with_site), main = 'No. of Cells Each Site is Observed In', breaks = 50)
hist(n_cells_with_site, main = 'No. of Cells Each Site is Observed In', breaks = 50)

sites_per_cell = colSums(f_binary_mat)
options(repr.plot.width = 8, repr.plot.height = 4)
par(mfrow = c(1, 2))
hist(log10(sites_per_cell), main = 'No. of Sites Observed per Cell', breaks = 50)
hist(sites_per_cell, main = 'No. of Sites Observed per Cell', breaks = 50)

pos = sapply(strsplit(rownames(f_binary_mat), split= ':'), tail , 1)
pos_len = sapply(strsplit(pos, split= '-'), function(x) as.numeric(x[2])-as.numeric(x[1]) )
par(mfrow = c(1, 1))
plot(pos_len, n_cells_with_site)
abline(v = f_binary_mat@Dim[2]*0.75)
dev.off()
f_binary_mat = f_binary_mat[ pos_len < 2000 & n_cells_with_site < f_binary_mat@Dim[2]*0.75, ]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels
genome(annotations) <- "mm10"

grange.counts <- StringToGRanges(rownames(f_binary_mat), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
f_binary_mat <- f_binary_mat[as.vector(grange.use),]

chrom_assay <- CreateChromatinAssay(
  counts = f_binary_mat,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = "/data/home/ruanlab/huangxingyu/scChAIR/ChIAV3/example/scATAC/ATAC10XlocalInOutPutATAC/fragment/8kmousecortex-1M.fragments.tsv.gz",
  min.cells = 1,
  min.features = 10,
  annotation = annotations
)