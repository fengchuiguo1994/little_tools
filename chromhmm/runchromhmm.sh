#samtools merge H3K4me3.bam ENCFF689TXP.bam ENCFF970WKG.bam
#samtools merge Control.bam ENCFF043PUQ.bam ENCFF948UNV.bam
#samtools merge H3K27me3.bam ENCFF630TXI.bam ENCFF692MXS.bam
#samtools merge H3K36me3.bam ENCFF129EML.bam ENCFF932DMH.bam
#samtools merge H3K4me1.bam ENCFF697MYF.bam ENCFF574QCS.bam
#samtools merge H3K27ac.bam ENCFF052BVG.bam ENCFF540LUD.bam
#samtools merge H3K9me3.bam ENCFF550RPT.bam ENCFF330FBC.bam
#samtools merge EP300.1.bam ENCFF848ZGM.bam ENCFF580PVH.bam
#samtools merge EP300.2.bam ENCFF359GMW.bam ENCFF336OIC.bam
#for i in EP300.2 EP300.1 H3K4me3 Control H3K27me3 H3K36me3 H3K4me1 H3K27ac H3K9me3
#do
#echo "picard MarkDuplicates I=$i.bam O=$i.rmd.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true M=$i.metrics REMOVE_DUPLICATES=true;samtools view -q 20 -o $i.flt.bam $i.rmd.bam;samtools index $i.flt.bam"|qsub -d . -N deal$i
#echo "bamCoverage -o $i.1x.bw --binSize 5 -b $i.flt.bam --numberOfProcessors 6 --minMappingQuality 15 --normalizeUsing RPGC --effectiveGenomeSize 2913022398"|qsub -d . -l nodes=1:ppn=6 -N bam2bw$i 

#done
 
#for i in H3K4me3 Control H3K27me3 H3K36me3 H3K4me1 H3K27ac H3K9me3
#do 
#macs2 callpeak -t ../A549.Pol2.DNA.bam -f BAM -g hs -n A549.Pol2.DNA --verbose 3 --trackline -B --SPMR --call-summits
#done


#java -jar ~/software/ChromHMM/ChromHMM/ChromHMM.jar BinarizeBam -b 200 ~/software/ChromHMM/ChromHMM/CHROMSIZES/hg38.txt ./ cellmarkfiletable.txt binary
#java -jar ~/software/ChromHMM/ChromHMM/ChromHMM.jar LearnModel binary 15stat/ 15 hg38


for i in EP300.1 EP300.2 H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3
do
bsub -q high -J run -n 1 -o macs.out -e macs.err -R span[hosts=1] "macs2 callpeak -t $i.flt.bam -c Control.flt.bam -f BAM -g 2913022398 -n $i --verbose 3 --trackline -B --SPMR"
bsub -q high -J run -n 1 -o macs.out -e macs.err -R span[hosts=1] "macs2 callpeak -t $i.flt.bam -c Control.flt.bam -f BAM -g 2913022398 -n $i --verbose 3 --trackline -B --SPMR --broad"
done
