library(grid)
library(RColorBrewer)

args=commandArgs(T)

mydat <- read.table(args[1],header=F)
# png(args[2], width = 800, height = 1000)
pdf(args[2], width = 8, height = 10)
plot_layout(24 * 3, 2, heights = unit(rep(c(0.75, 1.5, 0.75)/(24 * 3), 24), "npc"), widths = unit(c(0.1, 0.9), "npc"))
plot_chr_density(mydat)
dev.off()

## =================================================================== ##
## Set layout for plotting
## =================================================================== ##
plot_layout <- function(nrow, ncol, heights = NULL, widths = NULL) {
  if (is.null(heights)) {
    heights = unit(rep_len(1, nrow), "null")
  }
  if (is.null(widths)) {
    widths = unit(rep_len(1, ncol), "null")
  }
  grid.newpage()
  layout <- grid.layout(nrow = nrow, ncol = ncol,
                        heights = heights,
                        widths = widths)
  vp_overall <- viewport(name = "vp_overall", width = unit(0.95, "npc"), height = unit(0.95, "npc"), layout = layout)
  pushViewport(vp_overall)
}

## =================================================================== ##
## function to draw intra chromosomal interactions
## =================================================================== ##
plot_chr_density <- function(cyto) {
  cyto <- cyto[cyto[[1]] != "chrM", ]
  colnames(cyto) <- c("chrom", "start", "end", "value")
  # col_4 = colorRampPalette(brewer.pal(9, "Blues"))(18)     
  # col_4 = colorRampPalette(c("black", "white", "red"))(18)
  col_4 = colorRampPalette(c("Plum","lightblue","darkgreen","yellow", "orange", "red"))(18)
  
  basenumber=4
  cyto$col <- apply(cyto, 1, function(b){
    x = as.numeric(b[4])
    if(x <= 0){a = "white"}
    else if(x > 0*basenumber & x <=2*basenumber){a = col_4[1]}
    else if(x > 2*basenumber & x <=4*basenumber){a = col_4[2]}
    else if(x > 4*basenumber & x <=6*basenumber){a = col_4[3]}
    else if(x > 6*basenumber & x <=8*basenumber){a = col_4[4]}
    else if(x > 8*basenumber & x <=10*basenumber){a = col_4[5]}
    else if(x > 10*basenumber & x <=12*basenumber){a = col_4[6]}
    else if(x > 12*basenumber & x <=14*basenumber){a = col_4[7]}
    else if(x > 14*basenumber & x <=16*basenumber){a = col_4[8]}
    else if(x > 16*basenumber & x <=18*basenumber){a = col_4[9]}
    else if(x > 18*basenumber & x <=20*basenumber){a = col_4[10]}
    else if(x > 20*basenumber & x <=22*basenumber){a = col_4[11]}
    else if(x > 22*basenumber & x <=24*basenumber){a = col_4[12]}
    else if(x > 24*basenumber & x <=26*basenumber){a = col_4[13]}
    else if(x > 26*basenumber & x <=28*basenumber){a = col_4[14]}
    else if(x > 28*basenumber & x <=30*basenumber){a = col_4[15]}
    else if(x > 30*basenumber & x <=32*basenumber){a = col_4[16]}
    else if(x > 32*basenumber & x <=34*basenumber){a = col_4[17]}
    else{a = col_4[18]}})
    
  # write.table(cyto,file="test.out.txt",sep="\t",quote=F,row.names = F)
  # cyto$col <- "red"
  
  idx <- grepl("[0-9]+", cyto$chrom)
  chrom_no <- paste("chr", unique(c(sort(as.numeric(gsub("chr", "", cyto$chrom[idx]))), gsub("chr", "", sort(cyto$chrom[!idx])))), sep = "")
  cyto$chrom_no_f <- factor(cyto$chrom, levels = chrom_no)
  cyto <- cyto[with(cyto, order(chrom_no_f, start)), ]
  cyto <- cyto[, c("chrom", "start", "end", "col")]
  cyto[, c(1, 4)] <- sapply(cyto[, c(1, 4)], as.character)
  cyto[, 2:3] <- sapply(cyto[, 2:3], as.numeric)
  cyto_max <- max(cyto[["end"]])
  for (i in seq.int(length(chrom_no))) {
    if (grepl("X", chrom_no[i])) {
      chrom <- paste("chr", "X", sep = "")
    } else if (grepl("Y", chrom_no[i])) {
      chrom <- paste("chr", "Y", sep = "")
    }  else { chrom <- paste("chr", i, sep = "") }
    chrom_cyto <- cyto[which(cyto[[1]] == chrom), ]
    chrom_cyto <- chrom_cyto[order(chrom_cyto[[2]]), ]
    chrom_max <- max(chrom_cyto[["end"]])
    vp_chr = viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 1)
    grid.text(chrom, x = unit(0.5, "npc"), y = unit(0, "npc"),
              just = c("right", "bottom"), gp = gpar(col = "#323745", fontsize = 15) , vp = vp_chr)   
    vp = viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 2, name = "vp", xscale = c(0, cyto_max))
    chrom_cyto_2 <- rectGrob(x = unit(chrom_cyto[["start"]], "native"), y = unit(0, "npc"), width = unit((chrom_cyto[["end"]] - chrom_cyto[["start"]]) * 0.95, "native"), height = unit(1, "npc"),
                             just = c("left", "bottom"), gp = gpar(fill = chrom_cyto[["col"]], col = NA), vp = vp)           
    grid.draw(chrom_cyto_2)
    grid.rect(x = unit(0, "npc"), y = unit(1, "npc"), width = unit(chrom_max, "native"), height = unit(1, "npc"),
              just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)   

  }
}
