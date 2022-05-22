library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(pheatmap)
op <- options(stringsAsFactors=F)

mydata <- read.table("12356.txt",header=F,stringsAsFactors=F)
names(mydata)=c("chr","start","end","count","flag")
mydata.gr <- makeGRangesFromDataFrame(mydata)

myplot <- function(chrom="chr2",start=3700000,end=3800000){
    region <- data.frame(chr=c(chrom),start=c(start),end=c(end))
    region.gr <- makeGRangesFromDataFrame(region)
    x <- findOverlaps(region.gr,mydata.gr)

    res.gr <- mydata.gr[subjectHits(x),]
    res <- mydata[subjectHits(x),]

    step=ceiling((end-start)/100)
    bin=seq(from=start,to=end,by=step)
    binsize <- data.frame(chr=c(rep(chrom,100)),start=c(bin[1:100]),end=c(bin[2:101]),bin=c(seq(1:100)))
    binsize.gr <- makeGRangesFromDataFrame(binsize)
    x <- findOverlaps(binsize.gr,res.gr)

    aa=cbind(binsize[queryHits(x),],res[subjectHits(x),])
    mm=melt(aa,id=(c("flag","bin")))
    bb=dcast(mm,flag~bin)

    rownames(bb)=bb[,1]
    bb=bb[,-1]
    cluster=pheatmap(bb,silent=T,cluster_cols=FALSE,cluster_rows=TRUE)
    ord = rownames(bb[cluster$tree_row$order,])
    res$flag = factor(res$flag, levels=ord)

    ggplot(res) + 
    geom_segment(aes(x=start, xend=end, y=flag, yend=flag, color=flag), size=2) +
    geom_line(aes(x=start, y=flag,color=flag)) +
    theme_bw() +
    theme(legend.position="none") + 
    labs(title=paste(chrom,":",start,"-",end)) +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.y = element_blank()) + 
    theme(panel.grid =element_blank())
}
myplot(chrom="chr1",start=9400,end=10800)