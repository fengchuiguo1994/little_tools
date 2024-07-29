# python addBX.py test.R1.fastq.gz test.R2.fastq.gz SQ23016782
# pigz SQ23016782.BX.read1.fastq
# pigz SQ23016782.BX.read2.fastq
# fastqc -t 4 SQ23016782.BX.read1.fastq.gz SQ23016782.BX.read2.fastq.gz

i=SQ23016782
# sbatch -J trim.run -o trim.run.out -e trim.run.err -N 1 -n 10 --mem=10G --wrap=" java -jar /data/home/ruanlab/huangxingyu/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 10 ${i}.BX.read1.fastq.gz ${i}.BX.read2.fastq.gz ${i}.BX.read1.pair.fastq.gz ${i}.BX.read1.unpair.fastq.gz ${i}.BX.read2.pair.fastq.gz ${i}.BX.read2.unpair.fastq.gz ILLUMINACLIP:/data/home/ruanlab/huangxingyu/Tools/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:30 && fastqc -t 10 ${i}.BX.read1.pair.fastq.gz ${i}.BX.read1.unpair.fastq.gz ${i}.BX.read2.pair.fastq.gz ${i}.BX.read2.unpair.fastq.gz "

# mkdir align
# sbatch -J align.run -o align.run.out -e align.run.err -N 1 -n 10 --mem=10G --wrap=" bwa mem -t 10 /data/home/ruanlab/genome/hg38/hg38.fa ${i}.BX.read1.pair.fastq.gz ${i}.BX.read2.pair.fastq.gz | python renameTagPipe.py | samtools view -Sb -o align/${i}.BX.pair.bam - "
# sbatch -J align.run -o align.run.out -e align.run.err -N 1 -n 10 --mem=10G --wrap=" bwa mem -t 10 /data/home/ruanlab/genome/hg38/hg38.fa ${i}.BX.read1.unpair.fastq.gz | python renameTagPipe.py | samtools view -Sb -o align/${i}.BX.read1.unpair.bam - "
# sbatch -J align.run -o align.run.out -e align.run.err -N 1 -n 10 --mem=10G --wrap=" bwa mem -t 10 /data/home/ruanlab/genome/hg38/hg38.fa ${i}.BX.read2.unpair.fastq.gz | python renameTagPipe.py | samtools view -Sb -o align/${i}.BX.read2.unpair.bam - "

# samtools flagstat align/${i}.BX.pair.bam > align/${i}.BX.pair.bam.flagstat
# samtools flagstat align/${i}.BX.read1.unpair.bam > align/${i}.BX.read1.unpair.bam.flagstat
# samtools flagstat align/${i}.BX.read2.unpair.bam > align/${i}.BX.read2.unpair.bam.flagstat
# python bwaMAPQstat.py align/${i}.BX.pair.bam align/${i}.BX.pair.bam.MAPQ.stat
# python bwaMAPQstat.py align/${i}.BX.read1.unpair.bam align/${i}.BX.read1.unpair.bam.MAPQ.stat
# python bwaMAPQstat.py align/${i}.BX.read2.unpair.bam align/${i}.BX.read2.unpair.bam.MAPQ.stat

# sbatch -J align.run -o align.run.out -e align.run.err -N 1 -n 10 --mem=10G --wrap=" module load R-3.6.1 && java -cp CoLoci.jar RunBWA -i align/${i}.BX.pair.bam -r1 align/${i}.BX.read1.unpair.bam -r2 align/${i}.BX.read2.unpair.bam -o ${i} -m 30 -l 50 -e 500 -d 3000 -f 1 "

# zcat $i.all.frag.bed.gz | perl -lane '@tmp=split(/\t/,$_,3);@out=split(/;/,$tmp[2]);print join("\n",@out)' | sort --parallel=5 -k1,1 -k2,2n -k3,3n | bedtools genomecov -i - -bg -g /data/home/ruanlab/genome/hg38/hg38.fa.fai | sort --parallel=5 -k1,1 -k2,2n -k3,3n > $i.raw.bedgraph
# bedGraphToBigWig $i.raw.bedgraph /data/home/ruanlab/genome/hg38/hg38.fa.fai $i.raw.bw

# java -cp CoLoci.jar Mbed2HiC -i $i.all.frag.bed.gz -o $i.all.mbed2hic -j /data/home/ruanlab/huangxingyu/Tools/littletools/juicer_tools_1.22.01.jar -g /data/home/ruanlab/genome/hg38/hg38.fa.fai -s 1000 1>$i.all.mbed2hic.mbed2hic.out 2>$i.all.mbed2hic.err

# java -cp CoLoci.jar Mbed2HiC -i $i.samechrom.frag.bed.gz -o $i.samechrom.mbed2hic -j juicer_tools.1.8.9_jcuda.0.8.jar -g hg38.size -s 1000 1>$i.samechrom.mbed2hic.out 2>$i.samechrom.mbed2hic.err
# zcat SQ23016782.IGV.bed.gz | awk '$10>1' | awk -v OFS="\t" '{$9="255,0,0";print $0}' | sort --parallel=5 -k1,1 -k2,2n > SQ23016782.IGV.bed.sorted.bed12
