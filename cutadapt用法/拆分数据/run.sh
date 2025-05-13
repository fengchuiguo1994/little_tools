# cutadapt -j 10 -g file:linker.fa -n 50 -o test1.noLinker --info-file test1.Linker_info --discard -O 14 test1.fq
# cutadapt -j 10 -g file:linker.fa -n 50 -o test2.noLinker --info-file test2.Linker_info --discard -O 14 test2.fq

# cutadapt -j 10 -g file:linker.fa -n 50 -o testR1.noLinker --info-file testR1.Linker_info --discard -O 14 JZ25050738-Beam196-Beam196_combined_R1.fastq.gz
# cutadapt -j 10 -g file:linker.fa -n 50 -o testR2.noLinker --info-file testR2.Linker_info --discard -O 14 JZ25050738-Beam196-Beam196_combined_R2.fastq.gz

# python splitFile.py test1.Linker_info test2.Linker_info test1.fq test2.fq result
# python splitFilePro.py test1.Linker_info test2.Linker_info test1.fq test2.fq result
