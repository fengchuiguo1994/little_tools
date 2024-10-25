library(ClusterGVis)
mydat = read.table("C:\\Users\\dell\\Desktop\\表观修饰\\age.ATAC.txt", header=F)
mydat$peak = paste(mydat$V1, mydat$V2, mydat$V3, sep="-")
mydat = mydat[,c("peak", "V11", "V12", "V13")]
names(mydat) = c("peak","P95","P365","P730")
p95sum = sum(mydat$P95)
p365sum = sum(mydat$P365)
p730sum = sum(mydat$P730)

mydat$P95 = mydat$P95/p95sum * 1000000
mydat$P365 = mydat$P365/p365sum * 1000000
mydat$P730 = mydat$P730/p730sum * 1000000

exps = mydat

rownames(exps) = exps$peak
exps = exps[,-1]

# getClusters(exp = exps)
# cm <- clusterData(exp = exps, cluster.method = "mfuzz", cluster.num = 8)

# png("peak.cluster.png", width = 720)
# p1 = visCluster(object = cm, plot.type = "line")
# print(p1)
# dev.off()

# expsflt = exps[exps$P95 * 1.5 < exps$P365 | exps$P95 * 1.5  < exps$P730 | exps$P365 * 1.5  < exps$P730 | exps$P365 * 1.5 < exps$P95 | exps$P730 * 1.5  < exps$P95 | exps$P730 * 1.5  < exps$P365, ]
expsflt = exps[exps$P95 * 1.2 < exps$P365 | exps$P365 * 1.2  < exps$P730 | exps$P365 * 1.2 < exps$P95 | exps$P730 * 1.2  < exps$P365, ]
cm <- clusterData(exp = expsflt, cluster.method = "mfuzz", cluster.num = 16)
png("peak.cluster.flt.png", width = 720)
p1 = visCluster(object = cm, plot.type = "line")
print(p1)
dev.off()

up = expsflt[expsflt$P95 <= expsflt$P365 & expsflt$P365 <= expsflt$P730,]
cmup <- clusterData(exp = up, cluster.method = "mfuzz", cluster.num = 2)
png("peak.cluster.flt.up.png", width = 720)
p1 = visCluster(object = cmup, plot.type = "line")
print(p1)
dev.off()
write.table(cmup$wide.res, file="peak.up.peak.txt",quote=F,sep="\t")

down = expsflt[expsflt$P95 >= expsflt$P365 & expsflt$P365 >= expsflt$P730,]
cmdown <- clusterData(exp = down, cluster.method = "mfuzz", cluster.num = 2)
png("peak.cluster.flt.down.png", width = 720)
p1 = visCluster(object = cmdown, plot.type = "line")
print(p1)
dev.off()
write.table(cmdown$wide.res, file="peak.down.peak.txt",quote=F,sep="\t")

