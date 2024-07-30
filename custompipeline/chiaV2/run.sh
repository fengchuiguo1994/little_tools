# zcat ../SCG0192_GT22-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz | head -1000000 | gzip > SCG0192_GT22-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz
# zcat ../SCG0192_GT22-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz | head -1000000 | gzip > SCG0192_GT22-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz

bash chiaV2.sh chiaV2.config.sh SCG0192_GT22 . test
