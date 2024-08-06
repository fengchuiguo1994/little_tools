#cat ../02_mh_R1.paired.fq.gz ../03_mh_R1.paired.fq.gz > DNA_methy.1.fq.gz
#cat ../02_mh_R2.paired.fq.gz ../03_mh_R2.paired.fq.gz > DNA_methy.2.fq.gz
#echo "BatMeth2 pipel --fastp ~/miniconda2/bin/fastp -1 DNA_methy.1.fq.gz -2 DNA_methy.2.fq.gz -g ~/Genome/mh63S2/DNA_methy/MH63RS2.LNNK00000000.fasta -o combine -p 20 --gff=~/Genome/mh63S2/MH63RS2.LNNK00000000.v1.gff3 -f 1"|qsub -d . -l nodes=1:ppn=20 -N dnamathcombine

grep CHG combine.methratio.txt | awk '{print $1"\t"$2-1"\t"$2"\t"$7}' > CHG.bedgraph
grep CG combine.methratio.txt | awk '{print $1"\t"$2-1"\t"$2"\t"$7}' > CG.bedgraph
grep CHH combine.methratio.txt | awk '{print $1"\t"$2-1"\t"$2"\t"$7}' > CHH.bedgraph

bedGraphToBigWig CG.bedgraph ~/Allpet/xiaoqing/genome/MH63RS2.size CG.bw
bedGraphToBigWig CHG.bedgraph ~/Allpet/xiaoqing/genome/MH63RS2.size CHG.bw
bedGraphToBigWig CHH.bedgraph ~/Allpet/xiaoqing/genome/MH63RS2.size CHH.bw

computeMatrix scale-regions -S CG.bw CHG.bw CHH.bw -R ~/Allpet/xiaoqing/rloop/peak/all.peak -a 1000 -b 1000 -m 1000 --skipZeros --outFileName rloop.dnamethy.gz --numberOfProcessors 10


awk '{print $0"\t+"}' ~/Allpet/xiaoqing/rloop/peak/all.peak > rloop.bed
echo '/public/home/xyhuang/software/BatMeth2/bin/methyGff -o ./rloop -G /public/home/xyhuang/Genome/mh63S2/DNA_methy/MH63RS2.LNNK00000000.fasta -m combine.methratio.txt -b rloop.bed -B -P --TSS --TTS --GENE -s 0.025'|qsub -d .

cat ~/Allpet/xiaoqing/chrdseq/assemble/mh63.H3K4me3.cluster.FDRfiltered.txt | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3 | bedtools merge > chrd.bed
echo '/public/home/xyhuang/software/BatMeth2/bin/methyGff -o ./chrd -G /public/home/xyhuang/Genome/mh63S2/DNA_methy/MH63RS2.LNNK00000000.fasta -m combine.methratio.txt -b chrd.bed -B -P --TSS --TTS --GENE -s 0.025'|qsub -d .



echo '/public/home/xyhuang/software/BatMeth2/bin/methyGff -o ./none -G /public/home/xyhuang/Genome/mh63S2/DNA_methy/MH63RS2.LNNK00000000.fasta -m combine.methratio.txt -b ~/Allpet/xiaoqing/chrdseq/assemble/none.bed -B -P --TSS --TTS --GENE -s 0.025;Rscript ~/software/BatMeth2/scripts/methylevel.elements.r 0.025 none.Methylevel.1.txt none.test.pdf TSS TTS'|qsub -d .
echo '/public/home/xyhuang/software/BatMeth2/bin/methyGff -o ./chip -G /public/home/xyhuang/Genome/mh63S2/DNA_methy/MH63RS2.LNNK00000000.fasta -m combine.methratio.txt -b ~/Allpet/xiaoqing/chrdseq/assemble/chip.bed -B -P --TSS --TTS --GENE -s 0.025;Rscript ~/software/BatMeth2/scripts/methylevel.elements.r 0.025 chip.Methylevel.1.txt chip.test.pdf TSS TTS'|qsub -d .
echo '/public/home/xyhuang/software/BatMeth2/bin/methyGff -o ./rloop -G /public/home/xyhuang/Genome/mh63S2/DNA_methy/MH63RS2.LNNK00000000.fasta -m combine.methratio.txt -b ~/Allpet/xiaoqing/chrdseq/assemble/rloop.bed -B -P --TSS --TTS --GENE -s 0.025;Rscript ~/software/BatMeth2/scripts/methylevel.elements.r 0.025 rloop.Methylevel.1.txt rloop.test.pdf TSS TTS'|qsub -d .
echo '/public/home/xyhuang/software/BatMeth2/bin/methyGff -o ./both -G /public/home/xyhuang/Genome/mh63S2/DNA_methy/MH63RS2.LNNK00000000.fasta -m combine.methratio.txt -b ~/Allpet/xiaoqing/chrdseq/assemble/both.bed -B -P --TSS --TTS --GENE -s 0.025;Rscript ~/software/BatMeth2/scripts/methylevel.elements.r 0.025 both.Methylevel.1.txt both.test.pdf TSS TTS'|qsub -d .



echo 'computeMatrix scale-regions -S CG.bw CHG.bw CHH.bw -R /public/home/xyhuang/Allpet/xiaoqing/chrdseq/assemble/cis-trans.RNA.bed -a 2000 -b 2000 -m 2000 --skipZeros --outFileName cis-trans.RNA.bed.dnamethy.gz --numberOfProcessors 5 && plotProfile -m cis-trans.RNA.bed.dnamethy.gz --yMax 0.5 -out cis-trans.RNA.bed.dnamethy.gz.pdf --perGroup'|qsub -d . -l nodes=1:ppn=5
echo 'computeMatrix scale-regions -S CG.bw CHG.bw CHH.bw -R /public/home/xyhuang/Allpet/xiaoqing/chrdseq/assemble/trans.RNA.bed -a 2000 -b 2000 -m 2000 --skipZeros --outFileName trans.RNA.bed.dnamethy.gz --numberOfProcessors 5 && plotProfile -m trans.RNA.bed.dnamethy.gz --yMax 0.5 -out trans.RNA.bed.dnamethy.gz.pdf --perGroup'|qsub -d . -l nodes=1:ppn=5
echo 'computeMatrix scale-regions -S CG.bw CHG.bw CHH.bw -R /public/home/xyhuang/Allpet/xiaoqing/chrdseq/assemble/cis.RNA.bed -a 2000 -b 2000 -m 2000 --skipZeros --outFileName cis.RNA.bed.dnamethy.gz --numberOfProcessors 5 && plotProfile -m cis.RNA.bed.dnamethy.gz --yMax 0.5 -out cis.RNA.bed.dnamethy.gz.pdf --perGroup'|qsub -d . -l nodes=1:ppn=5

