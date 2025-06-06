# sed '1d' enhancer.txt | awk -F":" '{print $1"\t"$2}' | awk -F"-" '{print $1"\t"$2}' > enhancer.bed
# bedtools random -g ~/Tools/cellranger-7.1.0/genome/refdata-gex-mm10-2020-A/fasta/genome.fa.fai -l 1000 -n 50000 -seed 1234 | cut -f 1-3 > random.bed
# intersectBed -a ~/Haoxi20230215/finalresult/finalresult/1.peakcallingBrain/mouseBrain.raw.ATAC.final.bed -b enhancer.bed -v > other.ATAC.bed

# https://hgdownload.cse.ucsc.edu/goldenPath/mm10/phyloP60way/
# wget https://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw
# wget https://hgdownload.cse.ucsc.edu/goldenPath/mm10/phyloP4way/mm10.phyloP4way.bw

# conda activate deeptools
# computeMatrix scale-regions -S mm10.phastCons4way.bw --regionsFileName enhancer.bed other.ATAC.bed random.bed -a 3000 -b 3000 -m 3000 --skipZeros --outFileName mm10.phastCons4way.regions.gz -p 5
plotHeatmap -m mm10.phastCons4way.regions.gz -out mm10.phastCons4way.regions.Heatmap.pdf --perGroup
plotProfile -m mm10.phastCons4way.regions.gz -out mm10.phastCons4way.regions.Profile.pdf --perGroup
cp mm10.phastCons4way.regions.gz mm10.regions.1.gz
plotHeatmap -m mm10.regions.1.gz -out mm10.regions.1.Heatmap.pdf
plotProfile -m mm10.regions.1.gz -out mm10.regions.1.Profile.pdf

# computeMatrix scale-regions -S mm10.60way.phastCons.bw --regionsFileName enhancer.bed other.ATAC.bed random.bed -a 3000 -b 3000 -m 3000 --skipZeros --outFileName mm10.60way.phastCons.regions.gz -p 5
plotHeatmap -m mm10.60way.phastCons.regions.gz -out mm10.60way.phastCons.regions.Heatmap.pdf --perGroup
plotProfile -m mm10.60way.phastCons.regions.gz -out mm10.60way.phastCons.regions.Profile.pdf --perGroup
cp mm10.60way.phastCons.regions.gz mm10.regions.2.gz
plotHeatmap -m mm10.regions.2.gz -out mm10.regions.2.Heatmap.pdf
plotProfile -m mm10.regions.2.gz -out mm10.regions.2.Profile.pdf

# computeMatrix reference-point --referencePoint center -S mm10.phastCons4way.bw --regionsFileName enhancer.bed other.ATAC.bed random.bed -a 3000 -b 3000 --missingDataAsZero --outFileName mm10.phastCons4way.center.gz -p 10
plotHeatmap -m mm10.phastCons4way.center.gz -out mm10.phastCons4way.center.Heatmap.pdf --perGroup
plotProfile -m mm10.phastCons4way.center.gz -out mm10.phastCons4way.center.Profile.pdf --perGroup
cp mm10.phastCons4way.center.gz mm10.regions.3.gz
plotHeatmap -m mm10.regions.3.gz -out mm10.regions.3.Heatmap.pdf
plotProfile -m mm10.regions.3.gz -out mm10.regions.3.Profile.pdf

# computeMatrix reference-point --referencePoint center -S mm10.60way.phastCons.bw --regionsFileName enhancer.bed other.ATAC.bed random.bed -a 3000 -b 3000 --missingDataAsZero --outFileName mm10.60way.phastCons.center.gz -p 10
plotHeatmap -m mm10.60way.phastCons.center.gz -out mm10.60way.phastCons.center.Heatmap.pdf --perGroup
plotProfile -m mm10.60way.phastCons.center.gz -out mm10.60way.phastCons.center.Profile.pdf --perGroup
cp mm10.60way.phastCons.center.gz mm10.regions.4.gz
plotHeatmap -m mm10.regions.4.gz -out mm10.regions.4.Heatmap.pdf
plotProfile -m mm10.regions.4.gz -out mm10.regions.4.Profile.pdf