#!/bin/bash

NTHREAD=10 #cpu usage per job
ChIApipeV2="/data/home/ruanlab/huangxingyu/20230912/testtmp/chiaV2.sif"
# R1suf="_R1.fastq.gz" # suffix of the fastq files, modify if different
# R2suf="_R2.fastq.gz" # read2 suffix
R1suf="R1_001.fastq.gz"
R2suf="R3_001.fastq.gz"

### cpu parameters
genome="/data/home/ruanlab/huangxingyu/Haoxi20230315/20230205_pipeline/refs/mm10/bwa/mm10.fa"
linker=ACGTGATATTTCACGACTCT ## linker 3
minMapScore=30 #this is the default CPU memaln set
blacklist="null" # default: null

#cluster parameters:
addNone="T" # "T": PET from reads with linker and not. "F": PET from reads with linker.

#Chiapet params
extbp=500
selfbp=8000

#Peak calling params
genomelen='hs' #used in macs: human='hs' mouse='mm' fly='dm' or 2.7e9
macsq=0.001
peakext=150
shiftsize=$(( -$peakext/2 ))
allreads="F" # "T": will use all reads, "F": just use self lagation reads

#Chiasig params
PET=3
