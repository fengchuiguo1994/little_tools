# https://assets.nemoarchive.org/dat-1xp2hje

# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/QC.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/barcode.tsv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/cluster.annotation.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/cluster.membership.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/features.tsv.gz
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/matrix.mtx.gz
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/sample_metadata.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/umi_counts.h5

# mkdir 1xp2hje
# mv barcode.tsv features.tsv.gz matrix.mtx.gz 1xp2hje
# mv 1xp2hje/barcode.tsv 1xp2hje/barcodes.tsv
# gzip 1xp2hje/barcodes.tsv
# sed '1d' barcode.tsv | awk -F"," '{print $2}' | sed 's/\"//g' | gzip > barcodes.tsv.gz
