for i in Input1 Input3 S1 S2 S3 S4 S5 S6 SP9 SP10 SP11 SP12 SP13 SP14
do
    echo $i
    # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bowtie2.rmd.flt.bam -o ${i}.bowtie2.split -p dm3_ "
    # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bwa.rmd.flt.bam -o ${i}.bwa.split -p dm3_ "

    # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bowtie2.bam -o ${i}.bowtie2.raw.split -p dm3_ "
    # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bwa.bam -o ${i}.bwa.raw.split -p dm3_ "
done

for i in Input2 Input4 S7 S8 SP15 SP16
do
    echo $i
    # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bowtie2.rmd.flt.bam -o ${i}.bowtie2.split -p mm10_ "
    # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bwa.rmd.flt.bam -o ${i}.bwa.split -p mm10_ "

    # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bowtie2.bam -o ${i}.bowtie2.raw.split -p mm10_ "
    # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bwa.bam -o ${i}.bwa.raw.split -p mm10_ "
done


# cat *bowtie2.split.report.txt | grep -v n_exogenous | awk '{print $0"\t"1000000/$NF}' > all.bowtie2.split.report.txt
# cat *bwa.split.report.txt | grep -v n_exogenous | awk '{print $0"\t"1000000/$NF}' > all.bwa.split.report.txt
# spiker.py -t H3K27ac.sorted.bam -c control.sorted.bam --spikeIn --csf 1.23 --tsf 0.95 -o H3K27ac

# bsub -q q2680v2 -J ms -n 10 -o ms.out -e ms.err -R span[hosts=1] " multiBamSummary bins --bamfiles Input1.bowtie2.rmd.flt.bam Input3.bowtie2.rmd.flt.bam S1.bowtie2.rmd.flt.bam S2.bowtie2.rmd.flt.bam S3.bowtie2.rmd.flt.bam S4.bowtie2.rmd.flt.bam S5.bowtie2.rmd.flt.bam S6.bowtie2.rmd.flt.bam SP9.bowtie2.rmd.flt.bam SP10.bowtie2.rmd.flt.bam SP11.bowtie2.rmd.flt.bam SP12.bowtie2.rmd.flt.bam SP13.bowtie2.rmd.flt.bam SP14.bowtie2.rmd.flt.bam --minMappingQuality 20 --labels Input1 Input3 S1 S2 S3 S4 S5 S6 SP9 SP10 SP11 SP12 SP13 SP14 -out dm3.bowtie2.readCounts.npz --outRawCounts dm3.bowtie2.readCounts.tab --binSize 10000 --numberOfProcessors 10 "
# bsub -q q2680v2 -J ms -n 10 -o ms.out -e ms.err -R span[hosts=1] " multiBamSummary bins --bamfiles Input1.bwa.rmd.flt.bam Input3.bwa.rmd.flt.bam S1.bwa.rmd.flt.bam S2.bwa.rmd.flt.bam S3.bwa.rmd.flt.bam S4.bwa.rmd.flt.bam S5.bwa.rmd.flt.bam S6.bwa.rmd.flt.bam SP9.bwa.rmd.flt.bam SP10.bwa.rmd.flt.bam SP11.bwa.rmd.flt.bam SP12.bwa.rmd.flt.bam SP13.bwa.rmd.flt.bam SP14.bwa.rmd.flt.bam --minMappingQuality 20 --labels Input1 Input3 S1 S2 S3 S4 S5 S6 SP9 SP10 SP11 SP12 SP13 SP14 -out dm3.bwa.readCounts.npz --outRawCounts dm3.bwa.readCounts.tab --binSize 10000 --numberOfProcessors 10 "
# bsub -q q2680v2 -J ms -n 10 -o ms.out -e ms.err -R span[hosts=1] " multiBamSummary bins --bamfiles Input2.bowtie2.rmd.flt.bam Input4.bowtie2.rmd.flt.bam S7.bowtie2.rmd.flt.bam S8.bowtie2.rmd.flt.bam SP15.bowtie2.rmd.flt.bam SP16.bowtie2.rmd.flt.bam --minMappingQuality 20 --labels Input2 Input4 S7 S8 SP15 SP16 -out mm10.bowtie2.readCounts.npz --outRawCounts mm10.bowtie2.readCounts.tab --binSize 10000 --numberOfProcessors 10 "
# bsub -q q2680v2 -J ms -n 10 -o ms.out -e ms.err -R span[hosts=1] " multiBamSummary bins --bamfiles Input2.bwa.rmd.flt.bam Input4.bwa.rmd.flt.bam S7.bwa.rmd.flt.bam S8.bwa.rmd.flt.bam SP15.bwa.rmd.flt.bam SP16.bwa.rmd.flt.bam --minMappingQuality 20 --labels Input2 Input4 S7 S8 SP15 SP16 -out mm10.bwa.readCounts.npz --outRawCounts mm10.bwa.readCounts.tab --binSize 10000 --numberOfProcessors 10 "
# plotCorrelation -in dm3.bowtie2.readCounts.npz --corMethod spearman --plotTitle "Spearman" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o dm3.bowtie2.heatmap_SpearmanCorr_readCounts.pdf && plotCorrelation -in dm3.bowtie2.readCounts.npz --corMethod pearson --plotTitle "Pearson" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o dm3.bowtie2.heatmap_PearsonCorr_readCounts.pdf
# plotCorrelation -in mm10.bowtie2.readCounts.npz --corMethod spearman --plotTitle "Spearman" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o mm10.bowtie2.heatmap_SpearmanCorr_readCounts.pdf && plotCorrelation -in mm10.bowtie2.readCounts.npz --corMethod pearson --plotTitle "Pearson" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o mm10.bowtie2.heatmap_PearsonCorr_readCounts.pdf
# plotCorrelation -in dm3.bwa.readCounts.npz --corMethod spearman --plotTitle "Spearman" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o dm3.bwa.heatmap_SpearmanCorr_readCounts.pdf && plotCorrelation -in dm3.bwa.readCounts.npz --corMethod pearson --plotTitle "Pearson" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o dm3.bwa.heatmap_PearsonCorr_readCounts.pdf
# plotCorrelation -in mm10.bwa.readCounts.npz --corMethod spearman --plotTitle "Spearman" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o mm10.bwa.heatmap_SpearmanCorr_readCounts.pdf && plotCorrelation -in mm10.bwa.readCounts.npz --corMethod pearson --plotTitle "Pearson" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o mm10.bwa.heatmap_PearsonCorr_readCounts.pdf


# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K27ac.bowtie2.bam S1.bowtie2.rmd.flt.bam S2.bowtie2.rmd.flt.bam && samtools index SiHa.H3K27ac.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K4me1.bowtie2.bam S3.bowtie2.rmd.flt.bam S4.bowtie2.rmd.flt.bam  && samtools index SiHa.H3K4me1.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K4me3.bowtie2.bam S5.bowtie2.rmd.flt.bam S6.bowtie2.rmd.flt.bam && samtools index SiHa.H3K4me3.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.Sp1.bowtie2.bam S7.bowtie2.rmd.flt.bam S8.bowtie2.rmd.flt.bam && samtools index SiHa.Sp1.bowtie2.bam  "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K27ac.bowtie2.bam SP9.bowtie2.rmd.flt.bam SP10.bowtie2.rmd.flt.bam && samtools index SP.H3K27ac.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K4me1.bowtie2.bam SP11.bowtie2.rmd.flt.bam SP12.bowtie2.rmd.flt.bam && samtools index SP.H3K4me1.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K4me3.bowtie2.bam SP13.bowtie2.rmd.flt.bam SP14.bowtie2.rmd.flt.bam && samtools index SP.H3K4me3.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.Sp1.bowtie2.bam SP15.bowtie2.rmd.flt.bam SP16.bowtie2.rmd.flt.bam && samtools index SP.Sp1.bowtie2.bam "

# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K27ac.bwa.bam S1.bwa.rmd.flt.bam S2.bwa.rmd.flt.bam && samtools index SiHa.H3K27ac.bwa.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K4me1.bwa.bam S3.bwa.rmd.flt.bam S4.bwa.rmd.flt.bam  && samtools index SiHa.H3K4me1.bwa.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K4me3.bwa.bam S5.bwa.rmd.flt.bam S6.bwa.rmd.flt.bam && samtools index SiHa.H3K4me3.bwa.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.Sp1.bwa.bam S7.bwa.rmd.flt.bam S8.bwa.rmd.flt.bam && samtools index SiHa.Sp1.bwa.bam  "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K27ac.bwa.bam SP9.bwa.rmd.flt.bam SP10.bwa.rmd.flt.bam && samtools index SP.H3K27ac.bwa.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K4me1.bwa.bam SP11.bwa.rmd.flt.bam SP12.bwa.rmd.flt.bam && samtools index SP.H3K4me1.bwa.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K4me3.bwa.bam SP13.bwa.rmd.flt.bam SP14.bwa.rmd.flt.bam && samtools index SP.H3K4me3.bwa.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.Sp1.bwa.bam SP15.bwa.rmd.flt.bam SP16.bwa.rmd.flt.bam && samtools index SP.Sp1.bwa.bam "



for i in H3K27ac H3K4me1 H3K4me3
do
    for j in SiHa SP
    do
        echo $i
        # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${j}.${i}.bowtie2.bam -o ${j}.${i}.bowtie2.split -p dm3_ "
        # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${j}.${i}.bwa.bam -o ${j}.${i}.bwa.split -p dm3_ "
    done
done
for i in Sp1
do
    for j in SiHa SP
    do
        echo $i
        # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${j}.${i}.bowtie2.bam -o ${j}.${i}.bowtie2.split -p mm10_ "
        # bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${j}.${i}.bwa.bam -o ${j}.${i}.bwa.split -p mm10_ "
    done
done

# cat *.*.bowtie2.split.report.txt Input?.bowtie2.split.report.txt | grep -v n_exogenous | awk '{print $0"\t"1000000/$NF}' 
# SiHa.H3K27ac.bowtie2.bam	0	0	0	0	0	5915	53446975	3734602	0.267766
# SiHa.H3K4me1.bowtie2.bam	0	0	0	0	0	10170	45975680	6968612	0.143501
# SiHa.H3K4me3.bowtie2.bam	0	0	0	0	0	15537	50099972	10304583	0.0970442
# SiHa.Sp1.bowtie2.bam	0	0	0	0	0	31101	14845555	32470912	0.0307968
# SP.H3K27ac.bowtie2.bam	0	0	0	0	0	2374	58747591	1535134	0.651409
# SP.H3K4me1.bowtie2.bam	0	0	0	0	0	3325	52031762	2864862	0.349057
# SP.H3K4me3.bowtie2.bam	0	0	0	0	0	5087	59371899	4325836	0.231169
# SP.Sp1.bowtie2.bam	0	0	0	0	0	18642	43724447	11898767	0.0840423
# Input1.bowtie2.rmd.flt.bam	0	0	0	0	0	1914	9739531	1866272	0.535828
# Input2.bowtie2.rmd.flt.bam	0	0	0	0	0	6692	1915524	7028230	0.142283
# Input3.bowtie2.rmd.flt.bam	0	0	0	0	0	937	11306186	502054	1.99182
# Input4.bowtie2.rmd.flt.bam	0	0	0	0	0	7366	3867423	6079070	0.164499


# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K27ac.bowtie2.bam -c Input1.bowtie2.rmd.flt.bam --spikeIn --csf 0.535828 --tsf 0.267766 -o SiHa.H3K27ac.bowtie2.narrowpeak --bw --frip "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K4me1.bowtie2.bam -c Input1.bowtie2.rmd.flt.bam --spikeIn --csf 0.535828 --tsf 0.143501 -o SiHa.H3K4me1.bowtie2.narrowpeak --bw --frip "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K4me3.bowtie2.bam -c Input1.bowtie2.rmd.flt.bam --spikeIn --csf 0.535828 --tsf 0.0970442 -o SiHa.H3K4me3.bowtie2.narrowpeak --bw --frip "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.Sp1.bowtie2.bam -c Input2.bowtie2.rmd.flt.bam --spikeIn --csf 0.142283 --tsf 0.0307968 -o SiHa.Sp1.bowtie2.narrowpeak --bw --frip "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K27ac.bowtie2.bam -c Input3.bowtie2.rmd.flt.bam --spikeIn --csf 1.99182 --tsf 0.651409 -o SP.H3K27ac.bowtie2.narrowpeak --bw --frip "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K4me1.bowtie2.bam -c Input3.bowtie2.rmd.flt.bam --spikeIn --csf 1.99182 --tsf 0.349057 -o SP.H3K4me1.bowtie2.narrowpeak --bw --frip "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K4me3.bowtie2.bam -c Input3.bowtie2.rmd.flt.bam --spikeIn --csf 1.99182 --tsf 0.231169 -o SP.H3K4me3.bowtie2.narrowpeak --bw --frip "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.Sp1.bowtie2.bam -c Input4.bowtie2.rmd.flt.bam --spikeIn --csf 0.164499 --tsf 0.0840423 -o SP.Sp1.bowtie2.narrowpeak --bw --frip "




# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K27ac.bowtie2.bam -c Input1.bowtie2.rmd.flt.bam --spikeIn --csf 0.535828 --tsf 0.267766 -o SiHa.H3K27ac.bowtie2.broadpeak --bw --frip --broad "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K4me1.bowtie2.bam -c Input1.bowtie2.rmd.flt.bam --spikeIn --csf 0.535828 --tsf 0.143501 -o SiHa.H3K4me1.bowtie2.broadpeak --bw --frip --broad "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K4me3.bowtie2.bam -c Input1.bowtie2.rmd.flt.bam --spikeIn --csf 0.535828 --tsf 0.0970442 -o SiHa.H3K4me3.bowtie2.broadpeak --bw --frip --broad "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.Sp1.bowtie2.bam -c Input2.bowtie2.rmd.flt.bam --spikeIn --csf 0.142283 --tsf 0.0307968 -o SiHa.Sp1.bowtie2.broadpeak --bw --frip --broad "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K27ac.bowtie2.bam -c Input3.bowtie2.rmd.flt.bam --spikeIn --csf 1.99182 --tsf 0.651409 -o SP.H3K27ac.bowtie2.broadpeak --bw --frip --broad "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K4me1.bowtie2.bam -c Input3.bowtie2.rmd.flt.bam --spikeIn --csf 1.99182 --tsf 0.349057 -o SP.H3K4me1.bowtie2.broadpeak --bw --frip --broad "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K4me3.bowtie2.bam -c Input3.bowtie2.rmd.flt.bam --spikeIn --csf 1.99182 --tsf 0.231169 -o SP.H3K4me3.bowtie2.broadpeak --bw --frip --broad "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.Sp1.bowtie2.bam -c Input4.bowtie2.rmd.flt.bam --spikeIn --csf 0.164499 --tsf 0.0840423 -o SP.Sp1.bowtie2.broadpeak --bw --frip --broad "



for i in H3K27ac H3K4me1 H3K4me3 Sp1
do
    # bsub -q q2680v2 -J cal -n 10 -o cal.out -e cal.err -R span[hosts=1] " computeMatrix scale-regions -S SiHa.${i}.bowtie2.narrowpeak.treat.pileup.SpikeIn_scaled.bigWig SP.${i}.bowtie2.narrowpeak.treat.pileup.SpikeIn_scaled.bigWig -R hg19.HPV16.info -a 3000 -b 3000 -m 5000 --skipZeros --missingDataAsZero --outFileName hg19.HPV16.${i}.richheatmap.gz --numberOfProcessors 10 "
    plotHeatmap -m hg19.HPV16.${i}.richheatmap.gz -out hg19.HPV16.${i}.richheatmap.pdf --colorMap GnBu
    # bsub -q q2680v2 -J cal -n 10 -o cal.out -e cal.err -R span[hosts=1] " computeMatrix scale-regions -S SiHa.${i}.bowtie2.narrowpeak.treat.pileup.SpikeIn_scaled.bigWig SP.${i}.bowtie2.narrowpeak.treat.pileup.SpikeIn_scaled.bigWig -R ICG.gene.info -a 3000 -b 3000 -m 5000 --skipZeros --missingDataAsZero --outFileName ICG.gene.${i}.richheatmap.gz --numberOfProcessors 10 "
    plotHeatmap -m ICG.gene.${i}.richheatmap.gz -out ICG.gene.${i}.richheatmap.pdf --colorMap GnBu
done
