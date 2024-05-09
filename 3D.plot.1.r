library(plot3D)

mydat = read.table("img/combine.txt",header=F)
mydat$gex = log10((mydat[,5]+0.0001)/(mydat[,6]+0.001))
mydat$atac = log10((mydat[,12]+0.0001)/(mydat[,13]+0.001))
mydat$pet = log10((mydat[,19]+0.0001)/(mydat[,20]+0.001))
mydat[mydat$pet < -3,]$pet = -3
mydat$type = "mix"
mydat[mydat[,7]==mydat[,14] & mydat[,14]==mydat[,21] & mydat[,21]=="mm10",]$type = "mm10"
mydat[mydat[,7]==mydat[,14] & mydat[,14]==mydat[,21] & mydat[,21]=="hg38",]$type = "hg38"

tmp = mydat[,c("gex","atac","pet","type")]
colors = c("#E64B35FF","#D0DFE6FF","#4DBBD5FF")
colormap = c("#E64B35FF","#D0DFE6FF","#4DBBD5FF")
colors2 <- as.numeric(factor(tmp$type))

pdf("combine.txt.pdf",height=6,width=6, colormodel="cmyk")

with(tmp, scatter3D(gex, atac, pet, col=colormap,colvar = colors2, theta = 60, phi = 20,
                    xlab = "GEX", ylab = "ATAC", zlab = "PET",main = "hybrid",cex = 1,
                    bty = "f",ticktype = "detailed",d = 3,clab = c("type","class"),
                    adj = 0.5,font = 2,type = "s",pch = 16,box = TRUE,bg="black",alpha=0.4))

legend("right",title =  "omics",legend=c("hg38", "mix", "mm10"), pch=21,
       cex=1, y.intersp=2, pt.bg = colors, bg="white", bty="n")

dev.off()

