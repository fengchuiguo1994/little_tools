library(reshape2)
library(patchwork)
library(ggplot2)
args=commandArgs(T)
mydat = read.table(args[1]) # MOSCCA0063.raw.metadata.txt MOSCCA0063.raw.metadata.txt.pdf MOSCCA0063
mydat$cell = rownames(mydat)

is_decimal <- function(x) {
  # 判断是否为小数
  is_num <- is.numeric(x)
  is_int <- round(x) == x
  
  if (is_num && !is_int) {
    # 保留3位有效数字
    x <- format(round(x, 3), nsmall = 3)
    return(as.numeric(x))
  } else {
    return(x)
  }
}

aa = c('nCount_RNA','nFeature_RNA','percent.mt', 'percent.ribo','nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'PET_total', 'PET_inter', 'PET_intragt1000')
mydatflt = mydat[,c(aa[1], 'cell')]
mt = melt(mydatflt, id.vars = "cell")
med = is_decimal(median(mt$value))
mt$value = log10(mt$value)
bb = ggplot(mt, aes(x = variable , y = value)) +
       geom_violin() +
       stat_summary(fun = "median", geom = "text", label = med, position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
       theme_bw() +
       labs(x = args[3], y = "log10")

for (i in 2:length(aa)) {
  mydatflt = mydat[,c(aa[i], 'cell')]
  mt = melt(mydatflt, id.vars = "cell")
  med = is_decimal(median(mt$value))
  mt$value = log10(mt$value)
  cc = ggplot(mt, aes(x = variable , y = value)) +
    geom_violin() +
    stat_summary(fun = "median", geom = "text", label = med, position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
    theme_bw() +
    labs(x = args[3], y = "log10")
  bb = bb + cc
}
pdf(args[2])
print(bb)
dev.off()