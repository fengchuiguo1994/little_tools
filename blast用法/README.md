```
for i in {00..84}
do
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.$i.tar.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.$i.tar.gz.md5
    tar zxf nr.$i.tar.gz
done

wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar zxf taxdb.tar.gz

for i in {00..57}
do
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.$i.tar.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.$i.tar.gz.md5
    tar zxf nt.$i.tar.gz
done

```

```
samtools view -f 4 ~/20240422/hanjin.bam | awk '{print ">"$1"\n"$10}' | head -50000 > hanjin.unmap.fa
sbatch -p ruan_cpu -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 12 --mem=30G --wrap=" blastn -query hanjin.unmap.fa -out hanjin.unmap.output -db ../nt/nt -outfmt \"6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid scomname sblastname sskingdom staxids scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus ssciname sscinames\" -evalue 1e-5 -max_target_seqs 2 -num_threads 12 "


samtools view -f 4 -F 2048 OTNP4.bam | awk '{print ">"$1"\n"$10}' | head -50000 > OTNP4.unmap.fa
blastn -query OTTP4.unmap.fa -out OTTP4.unmap.output -db nt -outfmt "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid scomname sblastname sskingdom staxids scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus ssciname sscinames" -evalue 1e-5 -max_target_seqs 2 -num_threads 12
blastn -query hanjin.unmap.fa -out hanjin.unmap.output -db ../nt/nt -outfmt "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid scomname sblastname sskingdom staxids scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus ssciname sscinames" -evalue 1e-5 -max_target_seqs 2 -num_threads 12

perl getinfo.pl hanjin.unmap.output | cut -f 2 | sort | uniq -c | sed 's/^  *//;s/  */\t/' | sort -k1,1nr  > hanjin.unmap.tongji
```