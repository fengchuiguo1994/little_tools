# sbatch -p ruan_cpu -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 5 --mem=2G --wrap=" bamtofastq_linux --nthreads=5 --reads-per-fastq=500000000 possorted_genome_bam.bam test "

sbatch -p ruan_cpu -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 30 --mem=80G --wrap=" singularity run /data/home/ruanlab/huangxingyu/Haoxi20230315/20230205_pipeline/singularity/cellranger_6.0.0.sif count --transcriptome /data/home/ruanlab/huangxingyu/Tools/cellranger-7.1.0/genome/refdata-gex-mm10-2020-A --id 10xGEX --sample bamtofastq --fastqs=test/10xGEX_0_1_HJTGJDSX7 --localcores=30 --localmem=80 --include-introns --chemistry=ARC-v1 "

