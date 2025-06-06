# wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz
# wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz.md5
# wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi
# wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi.md5
# wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz
# wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz.md5
# wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi
# wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi.md5

# bcftools merge -0 -m snps -o patski.merge.vcf.gz C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz
# zcat patski.merge.vcf.gz | awk '{if(/^#/){print $0}else{if($7=="PASS"){print "chr"$0}}}' | gzip > patski.merge.rename.vcf.gz
# python parsevcf.py patski.merge.rename.vcf.gz patski.merge.hap.txt


# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 20 --mem=70G --wrap=" bwa mem -t 20 /data/home/ruanlab/huangxingyu/Haoxi20230315/20230205_pipeline/refs/mm10/bwa/mm10.fa /data/home/ruanlab/huangxingyu/Genome/phaseTest/testhicdata/patski.MboI.Rep3.R1_fastq.gz /data/home/ruanlab/huangxingyu/Genome/phaseTest/testhicdata/patski.MboI.Rep3.R2_fastq.gz | samtools view -@ 20 -Sb - | samtools sort -@ 20 -m 1500M -o patski.MboI.Rep3.bam - "
# sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 20 --mem=70G --wrap=" samtools index patski.MboI.Rep3.bam && bamCoverage -b patski.MboI.Rep3.bam -o patski.MboI.Rep3.bw -p 20 -bs 5 "
sbatch -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 20 --mem=70G --wrap=" samtools view -@ 20 -q 20 -o patski.MboI.Rep3.q20.bam patski.MboI.Rep3.bam && samtools index patski.MboI.Rep3.q20.bam "
