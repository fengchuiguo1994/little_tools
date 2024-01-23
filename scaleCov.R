mydat = read.table("k562.hg19.reptime.tsv.bedgraph", header = F)
mydat[mydat$V4 == -1,]$V4 = NA
mydat$V4 = scale(mydat$V4)
aaa = floor(min(mydat$V4,na.rm = T)) - 1
mydat[is.na(mydat$V4),]$V4 = aaa
write.table(mydat, file="k562.hg19.reptime.tsv.cg.bedgraph", sep="\t", quote=F, row.names = F, col.names = F)
