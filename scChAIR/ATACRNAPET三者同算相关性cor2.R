RNA=read.table("k562.RNA.txt",header=T)
ATAC=read.table("k562.ATAC.txt",header=T)
PET100K=read.table("k562.PET.100K.txt",header=T)
PETATAC=read.table("k562.PET.ATAC.txt",header=T)

rownames(RNA) = RNA[,1]
RNA = RNA[,-1]
rownames(ATAC) = ATAC[,1]
ATAC = ATAC[,-1]
rownames(PET100K) = PET100K[,1]
PET100K = PET100K[,-1]
rownames(PETATAC) = PETATAC[,1]
PETATAC = PETATAC[,-1]

RNA_gene <- rownames(RNA)
ATAC_gene <- rownames(ATAC)
PET100K_gene <- rownames(PET100K)
PETATAC_gene <- rownames(PETATAC)
gene <- intersect(intersect(intersect(RNA_gene,ATAC_gene),PET100K_gene),PETATAC_gene)
cor_RNAs <- RNA[gene,]
cor_ATACs <- ATAC[gene,]
cor_PET100K <- PET100K[gene,]
cor_PETATAC <- PETATAC[gene,]

cor_taylor <- function(X){
    tryCatch(
        { n <- ncol(X)
          return((1/sqrt(n))*sd(eigen(cor(X))$values)) },
        warning = function(w) {return(4)},
        error = function(e) { return(5)}
    )
}

cor_taylorspearman <- function(X){
    tryCatch(
        { n <- ncol(X)
          return((1/sqrt(n))*sd(eigen(cor(X,method="spearman"))$values)) },
        warning = function(w) {return(4)},
        error = function(e) { return(5)}
    )
}

cal = function(x) {
    hg = c(x)
    hg_RNA = cor_RNAs[hg,]
    hg_ATAC = cor_ATACs[hg,]
    hg_PET = cor_PETATAC[hg,]
    X=data.frame(RNA=t(hg_RNA),ATAC=t(hg_ATAC),PET=t(hg_PET))
    # print(head(X))
    a = cor_taylor(X)
    b = cor_taylorspearman(X)
    hg_PET = cor_PET100K[hg,]
    X=data.frame(RNA=t(hg_RNA),ATAC=t(hg_ATAC),PET=t(hg_PET))
    c = cor_taylor(X)
    d = cor_taylorspearman(X)
    print(sprintf("%s\t%f\t%f\t%f\t%f", x, a,b,c,d))
}

print("gene\tpromoterpearson\tpromoterspearman\t100Kpearson\t100Kspearman")
lapply(gene,cal)
