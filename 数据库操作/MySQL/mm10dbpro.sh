python uploadSQLPro.py addGenome -d mm10 -s scChAIR-viewer/mm10/mm10.size.txt
python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/mm10.chrband.cg.txt -d mm10 -p anno -n mm10.chrband -f chrband
awk -v OFS="\t" '{print $1,$2,$3,$4,$6,$7,$5,$8,"-"}' scChAIR-viewer/mm10/mm10.genome.max.bed > scChAIR-viewer/mm10/mm10.genome.max.v2.bed
python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/mm10.genome.max.v2.bed -d mm10 -p anno -n mm10.simple.anno -f anno
awk -v OFS="\t" '{print $1,$2,$3,$4,$6,$7,$5,$8,"-"}' scChAIR-viewer/mm10/mm10.genome.bed > scChAIR-viewer/mm10/mm10.genome.v2.bed
python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/mm10.genome.v2.bed -d mm10 -p anno -n mm10.full.anno -f anno

awk '{print $0"\tpeak"NR"\t0\t+\t-"}' scChAIR-viewer/mm10/patski.ATAC.G1.q001_peaks.final.bed > scChAIR-viewer/mm10/patski.ATAC.G1.q001_peaks.final.cg.v2.bed
python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/patski.ATAC.G1.q001_peaks.final.cg.v2.bed -d mm10 -p patski_scChAIR -n patski.ATAC.G1.peak -f bed6

python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/patski.ATAC.G1.RPKM.bedgraph -d mm10 -p patski_scChAIR -n patski.ATAC.G1.ATAC.RPKM.cov -f bedgraph

awk '{print $0"\t-"$4}' scChAIR-viewer/mm10/patski.GEX.G1.RPKM.bedgraph > scChAIR-viewer/mm10/patski.GEX.G1.RPKM.v2.ssbedgraph
python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/patski.GEX.G1.RPKM.v2.ssbedgraph -d mm10 -p patski_scChAIR -n patski.GEX.G1.RPKM.ssbedgraph -f ssbedgraph

python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/patski.PET.G1.clusters.gt2.bedpe -d mm10 -p patski_scChAIR -n patski.PET.G1.loop -f loop

python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/patski.GEX.G1.fragments.unsorted.tsv -d mm10 -p patski_scChAIR -n patski.scRNA.G1.scRNA -f scRNA

python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/patski.ATAC.G1.fragments.unsorted.tsv -d mm10 -p patski_scChAIR -n patski.scATAC.G1.scATAC -f scATAC

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,1}' scChAIR-viewer/mm10/patski.PET.bedpe.G1.bedpe > scChAIR-viewer/mm10/patski.PET.bedpe.G1.cg.v2.bedpe
python uploadSQLPro.py addTrack -i scChAIR-viewer/mm10/patski.PET.bedpe.G1.cg.v2.bedpe -d mm10 -p patski_scChAIR -n patski.scPET.G1.scPET -f scPET


awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,1}' scChAIR-viewer/mm10/patski.PET.bedpe.S.bedpe > scChAIR-viewer/mm10/patski.PET.bedpe.S.cg.v2.bedpe
python uploadSQLPro.py addMultiTrack -i scChAIR-viewer/mm10/patski.ATAC.S.RPKM.bedgraph,scChAIR-viewer/mm10/patski.GEX.S.RPKM.bedgraph,scChAIR-viewer/mm10/patski.PET.S.clusters.gt2.bedpe,scChAIR-viewer/mm10/patski.GEX.S.fragments.unsorted.tsv,scChAIR-viewer/mm10/patski.ATAC.S.fragments.unsorted.tsv,scChAIR-viewer/mm10/patski.PET.bedpe.S.cg.v2.bedpe -d mm10 -p patski_scChAIR -n patski_S -g cellcycleS -f ATACbedgraph,RNAbedgraph,loop,scRNA,scATAC,scPET

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,1}' scChAIR-viewer/mm10/patski.PET.bedpe.G2M.bedpe > scChAIR-viewer/mm10/patski.PET.bedpe.G2M.cg.v2.bedpe
python uploadSQLPro.py addMultiTrack -i scChAIR-viewer/mm10/patski.ATAC.G2M.RPKM.bedgraph,scChAIR-viewer/mm10/patski.GEX.G2M.RPKM.bedgraph,scChAIR-viewer/mm10/patski.PET.G2M.clusters.gt2.bedpe,scChAIR-viewer/mm10/patski.GEX.G2M.fragments.unsorted.tsv,scChAIR-viewer/mm10/patski.ATAC.G2M.fragments.unsorted.tsv,scChAIR-viewer/mm10/patski.PET.bedpe.G2M.cg.v2.bedpe -d mm10 -p patski_scChAIR -n patski_G2M -g cellcycleG2M -f ATACbedgraph,RNAbedgraph,loop,scRNA,scATAC,scPET


python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0088.sorted.4dn.hic.mcool -d mm10 -p patski_scChAIR -n SCG0088.mcool -f mcool -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0089.sorted.4dn.hic.mcool -d mm10 -p patski_scChAIR -n SCG0089.mcool -f mcool -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0090.sorted.4dn.hic.mcool -d mm10 -p patski_scChAIR -n SCG0090.mcool -f mcool -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0091.sorted.4dn.hic.mcool -d mm10 -p patski_scChAIR -n SCG0091.mcool -f mcool -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0092.sorted.4dn.hic.mcool -d mm10 -p patski_scChAIR -n SCG0092.mcool -f mcool -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0093.sorted.4dn.hic.mcool -d mm10 -p patski_scChAIR -n SCG0093.mcool -f mcool -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0088.sorted.4dn.hic -d mm10 -p patski_scChAIR -n SCG0088.hic -f hic -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0089.sorted.4dn.hic -d mm10 -p patski_scChAIR -n SCG0089.hic -f hic -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0090.sorted.4dn.hic -d mm10 -p patski_scChAIR -n SCG0090.hic -f hic -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0091.sorted.4dn.hic -d mm10 -p patski_scChAIR -n SCG0091.hic -f hic -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0092.sorted.4dn.hic -d mm10 -p patski_scChAIR -n SCG0092.hic -f hic -t 1
python uploadSQLPro.v2.py addTrack -i scChAIR-viewer/mm10/SCG0093.sorted.4dn.hic -d mm10 -p patski_scChAIR -n SCG0093.hic -f hic -t 1




# zcat scChAIR-viewer/hg38/mBC.bedpe.gz | awk '$1==$4' | awk -v OFS="\t" '{if($2<0){$2=0;} if($5<0){$5=0;} print $0}' | sed 's/,/_/' | awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7,1}' > scChAIR-viewer/hg38/mBC.cg.bedpe
# zcat scChAIR-viewer/hg38/hBM.bedpe.gz | awk '$1==$4' | awk -v OFS="\t" '{if($2<0){$2=0;} if($5<0){$5=0;} print $0}' | sed 's/,/_/' | awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7,1}' > scChAIR-viewer/hg38/hBM.cg.bedpe

# python splitbed.py scChAIR-viewer/hg38/hBM.cg.bedpe scChAIR-viewer/hg38/hBM.cg
# python splitbed.py scChAIR-viewer/hg38/mBC.cg.bedpe scChAIR-viewer/hg38/mBC.cg

# for i in {1..19} X; do python uploadSQLPro.v1.py addTrack -i scChAIR-viewer/hg38/mBC.cg.chr${i}.bedpe -d mm10 -p GAGE -n GAGE.mBC.scPET.chr${i} -f scPET; done
# for i in {1..19} X; do python uploadSQLPro.v1.py addTrack -i scChAIR-viewer/hg38/hBM.cg.chr${i}.bedpe -d hg38 -p GAGE -n GAGE.hBM.scPET.chr${i} -f scPET; done
