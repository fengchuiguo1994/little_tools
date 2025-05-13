args <- commandArgs(trailingOnly = TRUE)
library(pheatmap)
mydat = read.table(args[1],header=T)
# samples = c("CH001", "CH002", "CH005", "CH006", "CH008", "CH009", "CH010", "CH011", "CH012", "CH013", "CH014", "CH015", "CH016", "CH017", "CH018", "CH019", "CH020")
samples = c("CH001", "CH005", "CH008", "CH009", "CH011", "CH013", "CH015", "CH017", "CH019","CH002", "CH006", "CH010", "CH012", "CH014", "CH016", "CH018", "CH020")
rownames(mydat) = samples
colnames(mydat) = samples
pheatmap(mydat, cluster_rows = F, cluster_cols = F, number_format="%.2f", filename=args[2])
pheatmap(mydat, cluster_rows = T, cluster_cols = T, number_format="%.2f", filename=args[3])
# breaks = seq(0.6, 1, length.out = 101)
