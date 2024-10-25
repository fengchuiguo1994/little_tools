rm(list=ls())  #清空环境内变量
options(stringsAsFactors = F)  #避免自动将字符串转换为R语言因子
suppressMessages(library(GEOquery))
library(stringr)
library(ggplot2)
library(reshape2)
library(limma)


gset = getGEO('GSE240128', destdir=".", AnnotGPL = F, getGPL = F)
aa = gset$`GSE240128-GPL20795_series_matrix.txt.gz`
aa@experimentData
aa@phenoData@varMetadata
aa[["supplementary_file_1"]]
files = grep("supplementary_file", rownames(aa@phenoData@varMetadata), value=T)
bb = do.call(rbind, lapply(files, function(x) {aa[[x]]}))

cc = gset$`GSE240128-GPL21273_series_matrix.txt.gz`
cc@experimentData

cc@phenoData@varMetadata
cc[["supplementary_file_1"]]
files = grep("supplementary_file", rownames(cc@phenoData@varMetadata), value=T)
dd = do.call(rbind, lapply(files, function(x) {cc[[x]]}))

write.table(bb, file="bb.list", sep="\t", quote=F)
write.table(dd, file="dd.list", sep="\t", quote=F)
