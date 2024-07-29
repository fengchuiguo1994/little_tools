mydat = read.table("HCRS.results.flt.RawCounts",header=F)
names(mydat) = c("V1","V2","V3","HCRS0001","HCRS0002","HCRS0003")
cor(mydat[,c(4,5,6)])
cor(mydat[,c(4,5,6)],method="spearman")

mydat = read.table("HCRS.gene.count",header=T)
mydat = mydat[,c(7,8,9)]
names(mydat) = c("HCRS0001","HCRS0002","HCRS0003")
cor(mydat)
cor(mydat,method="spearman")

mydat = read.table("HCRS.fpkm",header=F)
mydat = mydat[,c(2,3,4)]
names(mydat) = c("HCRS0001","HCRS0002","HCRS0003")
cor(mydat)
cor(mydat,method="spearman")

mydat = read.table("HCRS.tpm",header=F)
mydat = mydat[,c(2,3,4)]
names(mydat) = c("HCRS0001","HCRS0002","HCRS0003")
cor(mydat)
cor(mydat,method="spearman")