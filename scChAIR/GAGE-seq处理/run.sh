# GSM7657698_contact_mBC-HiC-0716-1.pairs.gz
# GSM7657700_contact_mBC-HiC-0716-2.pairs.gz
# GSM7657702_contact_mBC-HiC-0814-1.pairs.gz
# cat <(zcat GSM7657698_contact_mBC-HiC-0716-1.pairs.gz | grep -v hg38 | sed 's/mm10_//g' | awk -v OFS="\t" '{print $1,$2-50,$2+50, $3,$4-50,$4+50, "mBC1_"$5, "read1"NR, 60, "+", "+"}') <(zcat GSM7657700_contact_mBC-HiC-0716-2.pairs.gz | grep -v hg38 | sed 's/mm10_//g' | awk -v OFS="\t" '{print $1,$2-50,$2+50, $3,$4-50,$4+50, "mBC2_"$5, "read2"NR, 60, "+", "+"}') <(zcat GSM7657702_contact_mBC-HiC-0814-1.pairs.gz | grep -v hg38 | sed 's/mm10_//g' | awk -v OFS="\t" '{print $1,$2-50,$2+50, $3,$4-50,$4+50, "mBC3_"$5, "read3"NR, 60, "+", "+"}') | pigz > mBC.bedpe



# GSM7657704_contact_CD-HiC-0722-1.pairs.gz
# GSM7657706_contact_CD-HiC-0722-2.pairs.gz
# GSM7657708_contact_CD-HiC-0722-4.pairs.gz
# cat <(zcat GSM7657704_contact_CD-HiC-0722-1.pairs.gz | grep -v mm10 | sed 's/hg38_//g' | awk -v OFS="\t" '{print $1,$2-50,$2+50, $3,$4-50,$4+50, "hBM1_"$5, "read1"NR, 60, "+", "+"}') <(zcat GSM7657706_contact_CD-HiC-0722-2.pairs.gz | grep -v mm10 | sed 's/hg38_//g' | awk -v OFS="\t" '{print $1,$2-50,$2+50, $3,$4-50,$4+50, "hBM2_"$5, "read2"NR, 60, "+", "+"}') <(zcat GSM7657708_contact_CD-HiC-0722-4.pairs.gz | grep -v mm10 | sed 's/hg38_//g' | awk -v OFS="\t" '{print $1,$2-50,$2+50, $3,$4-50,$4+50, "hBM3_"$5, "read3"NR, 60, "+", "+"}') | pigz > hBM.bedpe.gz


# grep hg38 K3-0208-PL2.info.txt | grep -v "well id" > GAGE.K562.list.txt
# zcat GSM7657692_contact_K3-0208-PL2-HiC.pairs.gz | awk '{if(NR==FNR){aa[$1]}else{if($5 in aa){print $0}}}' GAGE.K562.list.txt - | sed 's/hg38_//g' | awk -v OFS="\t" '{print $1,$2-50,$2+50, $3,$4-50,$4+50, $5, "read"NR, 60, "+", "+"}' > GAGE.K562.bedpe
# awk '{aa[$7]+=1}END{for(x in aa){print x"\t"aa[x]}}' GAGE.K562.bedpe > GAGE.K562.bedpe.stat


# awk '$9>=20' GAGE.K562.bedpe | awk -v OFS="\t" '{if($2<0){$2=0;} if($5<0){$5=0;} print $0}' | awk -v OFS="\t" '{print $1,$2,$3,$8,0,$10"\n"$4,$5,$6,$8,0,$11}' | bedtools bedtobam -i - -g ~/Genome/hg38/hg38.fa.fai > GAGE.K562.bam

# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=25G --wrap=" samtools sort -@ 10 -o GAGE.K562.sorted.bam GAGE.K562.bam && samtools index GAGE.K562.sorted.bam && bamCoverage -b GAGE.K562.sorted.bam -o GAGE.K562.bw -p 10 -bs 5 && bamCoverage -b GAGE.K562.sorted.bam -o GAGE.K562.RPKM.bw -p 10 -bs 5 --normalizeUsing RPKM "



# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=30G --wrap=" computeMatrix scale-regions -S GAGE.K562.bw --regionsFileName /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -a 2000 -b 2000 -m 2000 --missingDataAsZero --skipZeros --outFileName GAGE.K562.plotprofile.gz -p 10 "
# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=30G --wrap=" computeMatrix scale-regions -S GAGE.K562.RPKM.bw --regionsFileName /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -a 2000 -b 2000 -m 2000 --missingDataAsZero --skipZeros --outFileName GAGE.K562.plotprofile.RPKM.gz -p 10 "
# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=30G --wrap=" computeMatrix scale-regions -S GAGE.K562.bw --regionsFileName /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -a 2000 -b 2000 -m 2000 --missingDataAsZero --outFileName GAGE.K562.all.plotprofile.gz -p 10 "
# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=30G --wrap=" computeMatrix scale-regions -S GAGE.K562.RPKM.bw --regionsFileName /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -a 2000 -b 2000 -m 2000 --missingDataAsZero --outFileName GAGE.K562.all.plotprofile.RPKM.gz -p 10 "
# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=30G --wrap=" computeMatrix reference-point --referencePoint center -S GAGE.K562.bw --regionsFileName /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -a 5000 -b 5000 --missingDataAsZero --skipZeros --outFileName GAGE.K562.center.plotprofile.gz -p 10 "
# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=30G --wrap=" computeMatrix reference-point --referencePoint center -S GAGE.K562.RPKM.bw --regionsFileName /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -a 5000 -b 5000 --missingDataAsZero --skipZeros --outFileName GAGE.K562.center.plotprofile.RPKM.gz -p 10 "
# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=30G --wrap=" computeMatrix reference-point --referencePoint center -S GAGE.K562.bw --regionsFileName /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -a 5000 -b 5000 --missingDataAsZero --outFileName GAGE.K562.all.center.plotprofile.gz -p 10 "
# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=30G --wrap=" computeMatrix reference-point --referencePoint center -S GAGE.K562.RPKM.bw --regionsFileName /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -a 5000 -b 5000 --missingDataAsZero --outFileName GAGE.K562.all.center.plotprofile.RPKM.gz -p 10 "
ls *plotprofile*gz | while read i
do
    # plotHeatmap -m $i -out $i.heatmap.pdf --perGroup --zMin 0 --yMin 0
    # plotProfile -m $i -out $i.profile.pdf --perGroup --yMin 0
    sleep 0.001
done
# zcat GAGE.K562.all.center.plotprofile.gz | sed '1d' | awk -v sample="GAGE.K562" '{top=0;for(i=482;i<=531;i++){top+=$i;} left=0;for(i=7;i<=106;i++){left+=$i;} right=0;for(i=907;i<=1006;i++){right+=$i;} if(left+right==0){print "Inf\t"sample;} else{print 4*top/(left+right)"\t"sample}}' > GAGE.K562.all.center.plotprofile.gz.enrich.txt
# zcat GAGE.K562.all.center.plotprofile.RPKM.gz | sed '1d' | awk -v sample="GAGE.K562" '{top=0;for(i=482;i<=531;i++){top+=$i;} left=0;for(i=7;i<=106;i++){left+=$i;} right=0;for(i=907;i<=1006;i++){right+=$i;} if(left+right==0){print "Inf\t"sample;} else{print 4*top/(left+right)"\t"sample}}' > GAGE.K562.all.center.plotprofile.RPKM.gz.enrich.txt





# awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,"reads"NR,$7}' GAGE.K562.bedpe | grep -v -e chrM -e chrY > GAGE.K562.tmp.bedpe
# awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,"reads"NR,$7}' GAGE.K562.bedpe | grep -v -e chrM -e chrY | pairToBed -a - -b /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed -type both > GAGE.K562.PET.tmp.peak.both.bedpe
# awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,"reads"NR,$7}' GAGE.K562.bedpe | grep -v -e chrM -e chrY | pairToBed -a - -b /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/1.peakcalling/k562.flt.final.bed > GAGE.K562.PET.tmp.peak.bedpe
# pairToBed -a GAGE.K562.tmp.bedpe -b /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLine/hg38.promoterex3000.info > GAGE.K562.tmp.promoter.bedpe

# python /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/countTSSratio.py GAGE.K562.tmp.bedpe GAGE.K562.PET.tmp.peak.bedpe > GAGE.K562.PET.tmp.peak.bedpe.ratio
# python /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/countTSSratio.py GAGE.K562.tmp.bedpe GAGE.K562.PET.tmp.peak.both.bedpe > GAGE.K562.PET.tmp.peak.both.bedpe.ratio
# python /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/countTSSratio.py GAGE.K562.tmp.bedpe GAGE.K562.tmp.promoter.bedpe > GAGE.K562.tmp.promoter.bedpe.ratio
