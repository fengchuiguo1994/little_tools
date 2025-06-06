library(TreeDist)
library(pheatmap)
args=commandArgs(T)
pdf(args[2])
# mydat = read.table(args[1],header=T,check.names=F) # k562.raw.genome.txt
mydat = read.table('k562.raw.genome.txt',header=T,check.names=F) # k562.raw.genome.txt
rownames(mydat) = mydat[,1]
mydat = mydat[,-1]

normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

percentage <- function(x) {
  return (x / sum(x))
}

scalemydat = scale(mydat, center = TRUE, scale = TRUE) # Standard 标准化 # 中心化center：x-mean(x) # 标准化：中心化后除以标准差 (x-mean(x))/sd(x)
normalizemydat = sapply(mydat, normalize) # Normalize 归一化 0-1 (test <- preProcess(mydat, method = 'range')    testmydat <- predict(test, mydat)     table(testmydat == normalizemydat))
percentagemydat = sapply(mydat, percentage) # %
# scalemydat[1:15,1:15]

# (log2(50000) - 10) / 0.125 # bin44
# (log2(8000000) - 10) / 0.125 # bin104

# scalemydatflt = scalemydat[44:138,]
scalemydatflt = scalemydat[44:143,]
# normalizemydatflt = normalizemydat[44:138,]
normalizemydatflt = normalizemydat[44:143,]
# percentagemydatflt = percentagemydat[44:138,]
percentagemydatflt = percentagemydat[44:143,]

kl = KMeansPP(t(scalemydatflt),k=5,nstart=30)
# table(kl$cluster)
# spearmanmydat = cor(scalemydatflt, method="spearman")
# kl2 = KMeansPP(spearmanmydat,k=10,nstart=30)
# table(kl2$cluster)

saveRDS(kl,file="k562.raw.genome.txt.rds")
res = as.data.frame(kl$cluster)
names(res) = c("cluster")
resorder = res[order(res$cluster),]

tmydat = t(mydat)
tmydatorder = tmydat[rownames(resorder),]
mydatorder = t(tmydatorder)
pheatmap(mydatorder, cluster_rows=F, cluster_cols=F, scale="column", show_colnames = F)
write.table(res,file="k562.kmeans.result.txt",quote=F,sep="\t")
print(p1)

pdf("test.pdf",width=15)
# normalizemydatorder = sapply(mydatorder, normalize)
normalizemydatorder = normalizemydat[,rownames(resorder)]
pheatmap(normalizemydatorder, cluster_rows=F, cluster_cols=F, show_colnames = F)


pdf("test.pdf",width=15)
# percentagemydatorder = sapply(mydatorder, normalize)
percentagemydatorder = percentagemydat[,rownames(resorder)]
pheatmap(percentagemydatorder, cluster_rows=F, cluster_cols=F, show_colnames = F)

dev.off()
