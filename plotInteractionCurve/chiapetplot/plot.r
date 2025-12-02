#intra_overview

check_pkg <- function(pkg) {
  if(require(pkg, character.only = TRUE)){
    print(paste("Package", pkg, "is loaded correctly", sep = " "))
  } else {
    print(paste("Trying to install package", pkg, sep = " "))
    install.packages(pkg, repos="http://cran.us.r-project.org", dep = TRUE)
    if(require(pkg, character.only = TRUE)){
      print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
    } else{
      install.packages(pkg, repos="http://cran.rstudio.com/", dep = TRUE)
      if(require(pkg, character.only = TRUE)){
        print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
      } else{
        stop(paste("Couldn't install package", pkg, sep = " "));
      }
    }
  }
}
check_pkg("grid")
check_pkg("xtable")
check_pkg("RCircos")
## =================================================================== ##
## Plot intra chromosomal interaction
## =================================================================== ##
# h19 cytoband data downlaod
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
# mm10 cytoband data download
# http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/
#
setwd("~/科研/周强伟/项目/拟南芥-3D/chromcontact")
interactions <- read.table("h3k9me2.cluster.shuf4000.txt")
interactions_stat <- interactions[, 1:7]
interactions_stat[, c(1, 4)] <- sapply(interactions_stat[, c(1, 4)], as.character)
interactions_stat[, c(2, 3, 5, 6, 7)] <- sapply(interactions_stat[, c(2, 3, 5, 6, 7)], as.numeric)
interactions_same_chr <- interactions_stat[which(interactions_stat[[1]] == interactions_stat[[4]]), ]

#
## =================================================================== ##
## .PET_count_distribution.txt
## =================================================================== ##
pet_count <- read.table("ATCP-035_036_051_054-peak.PET_count_distribution.txt",sep = "\t", header = T)
#colnames(mapping_info) <- c("Title", "Information")
pet_count[[5]] <- round(pet_count[[5]], 4)
pet_count[[5]] <- paste(pet_count[[5]] * 100, "%", sep = "")
# get the PET conuts that have number of clusters more than 10000
pet_count[[2]] <- as.numeric(pet_count[[2]])
pet_10000 <- as.character(pet_count[-10, ][which(pet_count[[2]] > 20000), ][[1]])
#
source("Plotting_functions.R")

## =================================================================== ##
## Set color palettes for later use
## =================================================================== ##
fill_palette_1 <- c("#4D4D4D", "#C3C3C3")
fill_palette_2 <- c("#4D4D4D", "#969696", "#C3C3C3", "#E6E6E6", "#F5F5F5")
# fill_palette_3 <- c("#2A286B", "#4E69B2", "#EDA1AD", "#B9539F")
#fill_palette_3 <- c("#2A286B", "#4E69B2", "#72A7E5", "#E0ECF9")
fill_palette_3 <- c("#2F7ED8", "#72A7E5", "#95BDEB", "#DCE9F8")
fill_palette_4 <- c("#4D4D4D", "#969696", "#C3C3C3", "#E6E6E6")
format_axis_labels <- function(a) {
  sapply(a, function(x) {if (x >= 1000000) {x = paste(x/1000000, "M", sep = "")}
    else if (x>=1000) {x = paste(x/1000, "K", sep = "")}
    else {x = x}
  })
}
#
interactions_same_chr <- interactions_same_chr[which(!(interactions_same_chr[[7]] %in% pet_10000)), ]
clusters <- interactions_same_chr
peaks <- read.table("h3k9me2_peaks.shuf4000.bed")
cyto <- read.table("tair10_cytoBandIdeo.txt", fill = T, header = T)
peaks[[4]] <- (peaks[[4]] - min(peaks[[4]]))/(max(peaks[[4]]) - min(peaks[[4]]))

# filter cytoband data
cyto <- cyto[nchar(as.character(cyto[[1]])) <= 5, ]
#layout_nrow <- length(unique(cyto[[1]])) - 1
layout_nrow <- length(unique(cyto[[1]]))
#png("Rplot14.png", width = 800, height = 1000)
#pdf(paste(output_dir, "Rplot14.pdf", sep = ""), width = 10, height = 15)
#interactions_same_chr<- interactions_same_chr[interactions_same_chr$V7>4,]

pdf("h3k9me2.n4000.10M.plotchrcontact.pdf", width = 9, height = 11)
plot_layout(layout_nrow * 3, 2, heights = unit(rep(c(1.3, 0.4, 1.3)/(layout_nrow * 3), layout_nrow), "npc"), widths = unit(c(0.1, 0.9), "npc"))
plot_intra_chr_interaction(clusters, cyto, peaks)
dev.off()
