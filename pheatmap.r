args <- commandArgs(trailingOnly = TRUE)
library(pheatmap)
mydat = read.table(args[1],header=T)
pheatmap(mydat,display_numbers=T,number_format="%.2f",filename=args[2],breaks = seq(0.6, 1, length.out = 101))
