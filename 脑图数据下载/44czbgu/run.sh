# https://assets.nemoarchive.org/dat-ch1nqb7
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_Broad/QC.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_Broad/barcodes.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_Broad/cluster.annotation.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_Broad/cluster.membership.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_Broad/features.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_Broad/matrix.mtx.gz
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_Broad/sample_metadata.csv
# wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_Broad/umi_counts.h5

mkdir 44czbgu
mv matrix.mtx.gz 44czbgu
sed '1d' features.csv | awk '{print $1"\t"$1"\tGene Expression"}' | gzip > 44czbgu/features.tsv.gz
sed '1d' barcodes.csv | gzip > 44czbgu/barcodes.tsv.gz
