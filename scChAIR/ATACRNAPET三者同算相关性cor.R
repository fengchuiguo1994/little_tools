RNA=read.table("k562.RNA.txt",header=T)
ATAC=read.table("k562.ATAC.RNA.txt",header=T)
PET=read.table("k562.PET.txt",header=T)

rownames(RNA) = RNA[,1]
RNA = RNA[,-1]
rownames(ATAC) = ATAC[,1]
ATAC = ATAC[,-1]
rownames(PET) = PET[,1]
PET = PET[,-1]

RNA_gene <- rownames(RNA)
ATAC_gene <- rownames(ATAC)
gene <- intersect(RNA_gene,ATAC_gene)
cor_RNAs <- RNA[gene,]
cor_ATACs <- ATAC[gene,]
cor_PETs <- PET[gene,]

hg = c("MALAT1")
hg_RNA = cor_RNAs[hg,]
hg_ATAC = cor_ATACs[hg,]
hg_PET = PET[hg,]

cor_taylor <- function(X){
    # if(!inherits(X,"matrix")){
    #     stop("Input must be inherit ’matrix’ class.")
    # }
    tryCatch(
        { n <- ncol(X)
          return((1/sqrt(n))*sd(eigen(cor(X))$values)) },
        warning = function(w) {return(4)},
        error = function(e) { return(5)}
    )
}

tail(sort(rowSums(cor_RNAs)))

hg = c("MALAT1")
hg_RNA = cor_RNAs[hg,]
hg_ATAC = cor_ATACs[hg,]
hg_PET = PET[hg,]
X=data.frame(RNA=t(hg_RNA),ATAC=t(hg_ATAC),PET=t(hg_PET))
cor_taylor(X)


cal = function(x) {
    hg = c(x)
    hg_RNA = cor_RNAs[hg,]
    hg_ATAC = cor_ATACs[hg,]
    hg_PET = cor_PETs[hg,]
    X=data.frame(RNA=t(hg_RNA),ATAC=t(hg_ATAC),PET=t(hg_PET))
    return(cor_taylor(X))
}

cal("FTL")
lapply(gene,cal)

test = lapply(gene,cal)