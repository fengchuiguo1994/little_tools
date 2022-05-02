library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("gapit_functions.txt")
source("emma.txt")

#Step 1: Set working directory and import data
myY <- read.table("Japan.traits.txt", head = TRUE)
#myG <- read.table("mdp_genotype_test.hmp.txt" , head = FALSE)
myGD <- read.table("Japan2.numeric.txt", head = TRUE)
myGM <- read.table("Japan2.info.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT=GAPIT(Y=myY[,c(1,9)],GD=myGD,GM=myGM,PCA.total=0,group.from = 1,group.to = 1,group.by = 10,
  #sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  #sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  memo="ttest")

myGAPIT2=GAPIT(Y=myY[,c(1,9)],GD=myGD,GM=myGM,PCA.total=3,  group.from = 1,group.to = 1,group.by = 10,#sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  #sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  memo="GLM")

