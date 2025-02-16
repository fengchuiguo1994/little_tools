rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(ggplot2)
library(tidyverse)

# 单个基因查询
# GO通路
ego_ALL <- enrichGO(gene="EPCAM", OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="ALL", pvalueCutoff= 1,qvalueCutoff= 1, minGSSize = 1, maxGSSize = 100000)
ego_ALL@result
# KEGG通路
# 拿到这个基因的ENTREZID
bitr(gene="EPCAM", fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Hs.eg.db')
ekegg_ALL <- enrichKEGG(gene = "4072", organism = 'hsa', pvalueCutoff = 1, qvalueCutoff = 1,minGSSize = 1, maxGSSize = 100000)
ekegg_ALL@result

# MsigDB # https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp # msigdb.v2024.1.Hs.symbols.gmt
## 所有通路
geneset <- read.gmt("C:\\Users\\dell\\Desktop\\msigdb.v2024.1.Hs.symbols.gmt")
geneset_EPCAM <- geneset[geneset$gene == "EPCAM",]
head(geneset_EPCAM)
dim(geneset_EPCAM)
str(geneset_EPCAM)

