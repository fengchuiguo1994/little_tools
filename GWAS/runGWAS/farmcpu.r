library("bigmemory")
library("biganalytics")
library("compiler") 
source("gapit_functions.txt")
source("FarmCPU_functions.txt")

#Step 1: Set working directory and import data
myY <- read.table("../Japan.traits.txt", head = TRUE)
myGM <- read.table("../Japan2.info.txt", head = TRUE)
myGD <- read.big.matrix("../Japan2.numeric.txt", type="char", sep="\t", head = TRUE)
#Step 2: Run FarmCPU

myFarmCPU <- FarmCPU(
Y=myY[,c(1,9)],
GD=myGD,
GM=myGM
)

myP=as.numeric(myFarmCPU$GWAS$P.value)
myGI.MP=cbind(myGM[,-1],myP)
GAPIT.Manhattan(GI.MP=myGI.MP,seqQTN=mySim$QTN.position)
GAPIT.QQ(myP)

