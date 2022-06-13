args = commandArgs(T)
library(hicrep)
Nsample = length(args)
listsample = list()
x <- matrix(0,Nsample,Nsample)

for (i in 1:Nsample){
    listsample[[i]] = read.table(args[i],header=F)
}
for (i in 1:Nsample){
    for (j in 1:Nsample){
        tmp1 = listsample[[i]]
        tmp2 = listsample[[j]]
        h_hat <- htrain(tmp1, tmp2, 1000000, 5000000, 0:10)
        print(h_hat)
        processed <- prep(tmp1, tmp2, 1000000, h_hat, 5000000)
        scc.out <- get.scc(processed, 1000000, 5000000)
        x[i,j] = scc.out$scc[1]
    }
}

result = as.data.frame(x)
filename = basename(args)
colnames(result) = filename
rownames(result) = filename
write.table(result,file="result.out.txt",quote=F)

library(pheatmap)
bk <- c(seq(0.8,1,length=100))
pheatmap(result,filename="result.out.heatmap.pdf",breaks=bk,display_numbers=T,number_format="%.2f",number_color="black",fontsize_number=16)