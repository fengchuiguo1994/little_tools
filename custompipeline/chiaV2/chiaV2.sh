#!/usr/bin/bash

# all file must in your home.
configfile=$1 #configuration general to a set of libraries
fastqprefix=$2
outputdir=$3
outputprefix=$4

pairlabel='singlelinker.paired'
singlabel='singlelinker.single'
nonelabel='none'

source $configfile
scom="singularity run ${ChIApipeV2}"

mkdir -p ${outputdir}
LOGFILE=${outputdir}/${outputprefix}.log

####################### prepare input file
fastq1=${outputdir}/${outputprefix}_R1.fastq.gz
fastq2=${outputdir}/${outputprefix}_R2.fastq.gz
nf=$( ls ${fastqprefix}*${R1suf} | wc -l )
r1=$( ls ${fastqprefix}*${R1suf} )
r2=$( ls ${fastqprefix}*${R2suf} )

if [ "$nf" -gt 1 ]; then
   cat ${fastqprefix}*${R1suf} > $fastq1
   cat ${fastqprefix}*${R2suf} > $fastq2
else
    ln -s $r1 $fastq1
    ln -s $r2 $fastq2
fi

####################### perform linker detection and generation of different category of fastq files
echo "Linker detection " 1> ${LOGFILE}
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu stag -A ${linker} -W -T 18 -t ${NTHREAD} -O ${outputdir}/${outputprefix} $fastq1 $fastq2 1>> ${LOGFILE} 2>> ${LOGFILE}
echo "--- linker detection completed ----" 1>> ${LOGFILE}

####################### Get the stat
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu stat -s -p -T 18 -t 1 ${outputdir}/${outputprefix}.cpu 1>> ${outputdir}/${outputprefix}.stat 2>> ${LOGFILE}
echo "--- statistics done  ----" 1>> ${LOGFILE}

echo "--- pigziiping   ----" 1>> ${LOGFILE}
$scom pigz -p ${NTHREAD} ${outputdir}/${outputprefix}.singlelinker.paired.fastq 2>> ${LOGFILE}
$scom pigz -p ${NTHREAD} ${outputdir}/${outputprefix}.singlelinker.single.fastq 2>> ${LOGFILE}
$scom pigz -p ${NTHREAD} ${outputdir}/${outputprefix}.none.fastq 2>> ${LOGFILE}
$scom pigz -p ${NTHREAD} ${outputdir}/${outputprefix}.conflict.fastq 2>> ${LOGFILE}
$scom pigz -p ${NTHREAD} ${outputdir}/${outputprefix}.tied.fastq 2>> ${LOGFILE}
echo "--- pigziiping done  ----" 1>> ${LOGFILE}

####################### mapping
echo START ${outputdir}/${outputprefix} cpu memaln/pairing/span/deduplication/span .. 1>> ${LOGFILE}
for label in $pairlabel $singlabel $nonelabel
do
    echo "--- ${label} Mapping/pairing/deduplication ----" 1>> ${LOGFILE}
    echo "--- ${label} Mapping ----" 1>> ${LOGFILE}
    /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu memaln -T $minMapScore -t ${NTHREAD} ${genome} ${outputdir}/${outputprefix}.$label.fastq.gz 2>> ${LOGFILE} | $scom pigz -p ${NTHREAD} 1> ${outputdir}/${outputprefix}.$label.sam.gz
    $scom samtools view ${outputdir}/${outputprefix}.$label.sam.gz | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.$label.sam.gz.stat
    echo "--- ${label} Mapping done  ----" 1>> ${LOGFILE}

    echo "--- ${label} pairing ----" 1>> ${LOGFILE}
    /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu pair -s ${selfbp} -S -t ${NTHREAD} ${outputdir}/${outputprefix}.$label.sam.gz 1> ${outputdir}/${outputprefix}.$label.stat.xls 2>> ${LOGFILE}
    echo "--- ${label} pairing done  ----" 1>> ${LOGFILE}

    echo "--- ${label} deduplication ----" 1>> ${LOGFILE}
    if [ ${label} != ${singlabel} ];then
        if [ ${blacklist} != "null" ];then
            $scom samtools view -F 2048 -Sb ${outputdir}/${outputprefix}.$label.UxxU.bam | $scom bedtools intersect -nonamecheck -a - -b ${blacklist} -v | $scom samtools view -q ${minMapScore} -o ${outputdir}/${outputprefix}.$label.UxxU.rm.bam -
            $scom samtools view -F 2048 -Sb ${outputdir}/${outputprefix}.$label.UU.bam | $scom bedtools pairtobed -abam - -b ${blacklist} -type neither > ${outputdir}/${outputprefix}.$label.UU.rm.bam
            
        else
            $scom samtools view -F 2048 -Sb -o ${outputdir}/${outputprefix}.$label.UU.rm.bam ${outputdir}/${outputprefix}.$label.UU.bam 
            $scom samtools view -F 2048 -Sb -o ${outputdir}/${outputprefix}.$label.UxxU.rm.bam ${outputdir}/${outputprefix}.$label.UxxU.bam 
        fi
        $scom samtools view ${outputdir}/${outputprefix}.$label.UxxU.rm.bam | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.$label.UxxU.rm.bam.stat
        $scom samtools view ${outputdir}/${outputprefix}.$label.UU.rm.bam | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.$label.UU.rm.bam.stat
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu span -s ${selfbp} -g -t ${NTHREAD} ${outputdir}/${outputprefix}.$label.UU.rm.bam 2>> ${LOGFILE} 1> ${outputdir}/${outputprefix}.$label.UU.rm.span.xls

        # /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu dedup -g -t ${NTHREAD} ${outputdir}/${outputprefix}.$label.UU.rm.bam 1> ${outputdir}/${outputprefix}.$label.UU.rm.dedup.lc 2>> ${LOGFILE}
        # /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu dedup -g -t ${NTHREAD} ${outputdir}/${outputprefix}.$label.UxxU.rm.bam 1> ${outputdir}/${outputprefix}.$label.UxxU.rm.dedup.lc 2>> ${LOGFILE}
        $scom samtools sort -@ ${NTHREAD} -m 1G -o ${outputdir}/${outputprefix}.$label.UU.rm.csorted.bam ${outputdir}/${outputprefix}.$label.UU.rm.bam
        $scom samtools sort -@ ${NTHREAD} -m 1G -o ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.bam ${outputdir}/${outputprefix}.$label.UxxU.rm.bam
        $scom gatk MarkDuplicates --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M ${outputdir}/${outputprefix}.$label.UU.rm.csorted.MarkDuplicates.txt -O ${outputdir}/${outputprefix}.$label.UU.rm.csorted.nr.bam -I ${outputdir}/${outputprefix}.$label.UU.rm.csorted.bam 1>> ${LOGFILE} 2>&1
        $scom gatk MarkDuplicates --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.MarkDuplicates.txt -O ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.nr.bam -I ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.bam 1>> ${LOGFILE} 2>&1
        $scom samtools sort -n -@ ${NTHREAD} -o ${outputdir}/${outputprefix}.$label.UU.rm.nr.bam ${outputdir}/${outputprefix}.$label.UU.rm.csorted.nr.bam
        $scom samtools sort -n -@ ${NTHREAD} -o ${outputdir}/${outputprefix}.$label.UxxU.rm.nr.bam ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.nr.bam

        $scom samtools view ${outputdir}/${outputprefix}.$label.UU.rm.nr.bam | awk -F":" '{print $1}' | uniq >  ${outputdir}/${outputprefix}.$label.UU.rm.nr.bam.stat
        $scom samtools view ${outputdir}/${outputprefix}.$label.UxxU.rm.nr.bam | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.$label.UxxU.rm.nr.bam.stat

        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu span -s ${selfbp} -g -t ${NTHREAD} ${outputdir}/${outputprefix}.$label.UU.rm.nr.bam 2>> ${LOGFILE} 1> ${outputdir}/${outputprefix}.$label.UU.rm.nr.span.xls

        if [ ${label} == ${nonelabel} ];then
            $scom filterinsertsizeless1000 ${outputdir}/${outputprefix}.$label.UU.rm.nr.bam ${outputdir}/${outputprefix}.$label.UU.rm.nr
        fi
    else
        if [ ${blacklist} != "null" ];then
            $scom samtools view -F 2048 -Sb ${outputdir}/${outputprefix}.$label.UxxU.bam | $scom bedtools intersect -nonamecheck -a - -b ${blacklist} -v | $scom samtools view -q ${minMapScore} -o ${outputdir}/${outputprefix}.$label.UxxU.rm.bam -
        else
            $scom samtools view -F 2048 -Sb -o ${outputdir}/${outputprefix}.$label.UxxU.rm.bam ${outputdir}/${outputprefix}.$label.UxxU.bam
        fi
        $scom samtools view ${outputdir}/${outputprefix}.$label.UxxU.rm.bam | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.$label.UxxU.rm.bam.stat

        # /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu dedup -g -t ${NTHREAD} ${outputdir}/${outputprefix}.$label.UxxU.rm.bam 1> ${outputdir}/${outputprefix}.$label.UxxU.rm.dedup.lc 2>> ${LOGFILE}
        $scom samtools sort -@ ${NTHREAD} -m 1G -o ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.bam ${outputdir}/${outputprefix}.$label.UxxU.rm.bam
        $scom gatk MarkDuplicates --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.MarkDuplicates.txt -O ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.nr.bam -I ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.bam 1>> ${LOGFILE} 2>&1
        $scom samtools sort -n -@ ${NTHREAD} -o ${outputdir}/${outputprefix}.$label.UxxU.rm.nr.bam ${outputdir}/${outputprefix}.$label.UxxU.rm.csorted.nr.bam
        $scom samtools view ${outputdir}/${outputprefix}.$label.UxxU.rm.nr.bam | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.$label.UxxU.rm.nr.bam.stat
    fi
    echo "--- ${label} deduplication done  ----" 1>> ${LOGFILE}
done

####################### peak calling
echo "--- peak calling ----" 1>> ${LOGFILE}
$scom samtools merge -@ ${NTHREAD} -o ${outputdir}/${outputprefix}.csorted.bam ${outputdir}/${outputprefix}.$pairlabel.UU.rm.csorted.nr.bam ${outputdir}/${outputprefix}.$pairlabel.UxxU.rm.csorted.nr.bam ${outputdir}/${outputprefix}.$nonelabel.UU.rm.csorted.nr.bam ${outputdir}/${outputprefix}.$nonelabel.UxxU.rm.csorted.nr.bam ${outputdir}/${outputprefix}.$singlabel.UxxU.rm.csorted.nr.bam
samtools index ${outputdir}/${outputprefix}.csorted.bam
if [ ${allreads} == "F" ];then
    $scom samtools sort -@ ${NTHREAD} -m 1G -o ${outputdir}/${outputprefix}.formacs2.bam ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.singleton.bam
else
    ln -s ${outputdir}/${outputprefix}.csorted.bam ${outputdir}/${outputprefix}.formacs2.bam
fi
$scom macs2 callpeak --keep-dup all --nomodel --shift $shiftsize  --extsize $peakext -B --SPMR -t ${outputdir}/${outputprefix}.formacs2.bam -g $genomelen -n ${outputdir}/${outputprefix} --qvalue $macsq >> ${LOGFILE} 2>> ${LOGFILE}
echo "--- peak calling done  ----" 1>> ${LOGFILE}

####################### cluster
echo "--- cluster ----" 1>> ${LOGFILE}
$scom samtools merge -n -@ ${NTHREAD} -o ${outputdir}/${outputprefix}.pet.bam ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.bam ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.pet.bam
$scom samtools view ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.bam | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.bam.stat
$scom samtools view ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.pet.bam | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.pet.bam.stat

$scom cpu span -s ${selfbp} -g -t ${NTHREAD} ${outputdir}/${outputprefix}.pet.bam 2>> ${LOGFILE} 1> ${outputdir}/${outputprefix}.pet.span.xls
$scom cpu span -s ${selfbp} -g -t ${NTHREAD} ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.bam 2>> ${LOGFILE} 1> ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.span.xls
$scom cpu span -s ${selfbp} -g -t ${NTHREAD} ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.pet.bam 2>> ${LOGFILE} 1> ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.pet.span.xls

$scom samtools view ${outputdir}/${outputprefix}.pet.bam | awk -F":" '{print $1}' | uniq > ${outputdir}/${outputprefix}.pet.bam.stat
if [ ${addNone} == "T" ];then
    ln -s ${outputdir}/${outputprefix}.pet.bam ${outputdir}/${outputprefix}.forcluster.bam
else
    ln -s ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.bam ${outputdir}/${outputprefix}.forcluster.bam
fi
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom cpu cluster -s $selfbp -M -B 1000 -5 5,-20 -3 5,480 -t ${NTHREAD} -O ${outputdir}/${outputprefix}.e$extbp -j -x -v 1 -g ${outputdir}/${outputprefix}.forcluster.bam 2>> ${LOGFILE} 1>> ${LOGFILE}
zcat ${outputdir}/${outputprefix}.e$extbp.clusters.cis.chiasig.gz | awk '$7+0.0 > 1' > ${outputdir}/${outputprefix}.e$extbp.clusters.cis.chiasig.BE2
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M' $scom ChiaSig -c $PET -t $NTHREAD -p -m $selfbp ${outputdir}/${outputprefix}.e$extbp.clusters.cis.chiasig.BE2 >> ${LOGFILE} 2>> ${LOGFILE}

$scom bedtools bamtobed -i ${outputdir}/${outputprefix}.forcluster.bam -bedpe > ${outputdir}/${outputprefix}.forcluster.bedpe
$scom bedpe2ipet ${outputdir}/${outputprefix}.forcluster.bedpe ${outputdir}/${outputprefix}.forcluster.ipet
$scom loopCallingGivenAnchor ${outputdir}/${outputprefix}_peaks.narrowPeak ${outputdir}/${outputprefix}.forcluster.ipet ${outputdir}/${outputprefix}.forcluster.ipet.loops > ${outputdir}/${outputprefix}.forcluster.ipet.loops.log 2>&1
echo "--- cluster done  ----" 1>> ${LOGFILE}

####################### for viewer
$scom samtools view -H ${outputdir}/${outputprefix}.csorted.bam | grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? > ${outputdir}/${outputprefix}.genome.length
$scom bedClip ${outputdir}/${outputprefix}_treat_pileup.bdg ${outputdir}/${outputprefix}.genome.length ${outputdir}/${outputprefix}_treat_pileup.clip.bdg
sort --parallel=${NTHREAD} -k1,1 -k2,2n ${outputdir}/${outputprefix}_treat_pileup.clip.bdg > ${outputdir}/${outputprefix}_treat_pileup.clip.sorted.bdg
$scom bedGraphToBigWig ${outputdir}/${outputprefix}_treat_pileup.clip.sorted.bdg ${outputdir}/${outputprefix}.genome.length ${outputdir}/${outputprefix}.treat_pileup.NDP.bw

$scom bamCoverage -bs 5 -b ${outputdir}/${outputprefix}.csorted.bam -o ${outputdir}/${outputprefix}.forBASIC.bw -p ${NTHREAD} 1>> ${LOGFILE} 2>&1
$scom bamCoverage -bs 5 -b ${outputdir}/${outputprefix}.csorted.bam -o ${outputdir}/${outputprefix}.forBASIC.RPKM.bw -p ${NTHREAD} --normalizeUsing RPKM 1>> ${LOGFILE} 2>&1
$scom bamCoverage -bs 5 -b ${outputdir}/${outputprefix}.csorted.bam -o ${outputdir}/${outputprefix}.forBASIC.bedgraph -p ${NTHREAD} -of bedgraph 1>> ${LOGFILE} 2>&1
$scom bamCoverage -bs 5 -b ${outputdir}/${outputprefix}.csorted.bam -o ${outputdir}/${outputprefix}.forBASIC.RPKM.bedgraph -p ${NTHREAD} --normalizeUsing RPKM  -of bedgraph 1>> ${LOGFILE} 2>&1

$scom juicer_tools pre -r 2500000,2000000,1000000,500000,250000,200000,100000,50000,25000,20000,10000,5000,2500,2000,1000 ${outputdir}/${outputprefix}.e${extbp}.juice.gz ${outputdir}/${outputprefix}.hic ${outputdir}/${outputprefix}.genome.length 2>> ${LOGFILE} 1>> ${LOGFILE}

$scom bedpe2ipetle1000 ${outputdir}/${outputprefix}.forcluster.bedpe ${outputdir}/${outputprefix}.forcluster.gt1000.ipet
cat ${outputdir}/${outputprefix}.forcluster.gt1000.ipet | awk -v OFS="\t" '{print 0,$1,int(($2+$3)/2),0,0,$4,int(($5+$6)/2),1}' | sort --parallel=${NTHREAD} -k2,2 -k6,6 -k3,3n -k7,7n > ${outputdir}/${outputprefix}.forcluster.gt1000.ipet.4dn 
$scom juicer_tools pre -r 2500000,2000000,1000000,500000,250000,200000,100000,50000,25000,20000,10000,5000,2500,2000,1000 ${outputdir}/${outputprefix}.forcluster.gt1000.ipet.4dn ${outputdir}/${outputprefix}.forcluster.gt1000.ipet.4dn.hic ${outputdir}/${outputprefix}.genome.length 2>> ${LOGFILE} 1>> ${LOGFILE}

#Get the summary
out_file=${outputdir}/${outputprefix}.final_stats.tsv
n_pet=$( cat ${outputdir}/${outputprefix}.stat | grep "Total pairs" | awk -F'[ \t]' '{print $3}' )
echo -e "Total_PE_reads\t"${n_pet} > ${out_file}
read_link=$( cat ${outputdir}/${outputprefix}.stat | grep "Linker detected" | awk -F '[ \t]' '{print $3}' | xargs printf "%'.f")
echo -e "Linker_detected_in_pair_reads\t"${read_link} >> ${out_file}
pet_link=$( cat ${outputdir}/${outputprefix}.stat | grep "Single Linker 2 tags" | awk -F '[ \t]' '{print $6}' | xargs printf "%'.f")
echo -e "PET_with_linker\t"${pet_link} >> ${out_file}

echo "" >> ${out_file}
pairlabeltotal=$(cat ${outputdir}/${outputprefix}.$pairlabel.sam.gz.stat | wc -l | xargs printf "%'.f")
echo -e "paired_reads\t"${pairlabeltotal} >> ${out_file}
pairUUuniq=$(cat ${outputdir}/${outputprefix}.$pairlabel.UU.rm.bam.stat | wc -l | xargs printf "%'.f")
pairUxxUuniq=$(cat ${outputdir}/${outputprefix}.$pairlabel.UxxU.rm.bam.stat | wc -l | xargs printf "%'.f")
echo -e "paired_UU_reads\t"${pairUUuniq} >> ${out_file}
echo -e "paired_UxxU_reads\t"${pairUxxUuniq} >> ${out_file}
pairUUuniq=$(cat ${outputdir}/${outputprefix}.$pairlabel.UU.rm.bam.stat | wc -l)
pairUxxUuniq=$(cat ${outputdir}/${outputprefix}.$pairlabel.UxxU.rm.bam.stat | wc -l)
pairUUuniqnr=$(cat ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.bam.stat | wc -l | xargs printf "%'.f")
pairUxxUuniqnr=$(cat ${outputdir}/${outputprefix}.$pairlabel.UxxU.rm.nr.bam.stat | wc -l | xargs printf "%'.f")
echo -e "paired_UU_reads_dedup\t"${pairUUuniqnr} >> ${out_file}
echo -e "paired_UxxU_reads_dedup\t"${pairUxxUuniqnr} >> ${out_file}
pairUUuniqnr=$(cat ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.bam.stat | wc -l )
pairUxxUuniqnr=$(cat ${outputdir}/${outputprefix}.$pairlabel.UxxU.rm.nr.bam.stat | wc -l )
pair_self_lig=$( cat ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.span.xls | grep "second/best<0.95" -A5 |  awk -F '[\t]' '{if(NR==4)print $2}' )
pair_self_lig=$( printf "%'.f" ${pair_self_lig} )
echo -e "linker_Self-ligation_PET\t"${pair_self_lig} >> ${out_file}
pair_intra_chr_pet=$( cat ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==5)print $2}' )
pair_inter_chr_pet=$( cat ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==2)print $2}' )
pair_inter_lig_all=$( echo "${pair_intra_chr_pet} + ${pair_inter_chr_pet}" | bc )
pair_inter_lig_all=$( printf "%'.f" ${pair_inter_lig_all} )
echo -e "linker_Inter-ligation_PET\t"${pair_inter_lig_all} >> ${out_file}

pair_intra_chr_pet=$( printf "%'.f" ${pair_intra_chr_pet} )
echo -e "linker_Intra-chr_PET\t"${pair_intra_chr_pet} >> ${out_file}
pair_inter_chr_pet=$( printf "%'.f" ${pair_inter_chr_pet} )
echo -e "linker_Inter-chr_PET\t"${pair_inter_chr_pet} >> ${out_file}

echo "" >> ${out_file}
singlabeltotal=$(cat ${outputdir}/${outputprefix}.$singlabel.sam.gz.stat | wc -l | xargs printf "%'.f")
echo -e "single_reads\t"${singlabeltotal} >> ${out_file}
singUxxUuniq=$(cat ${outputdir}/${outputprefix}.$singlabel.UxxU.rm.bam.stat | wc -l | xargs printf "%'.f")
echo -e "single_UxxU_reads\t"${singUxxUuniq} >> ${out_file}
singUxxUuniq=$(cat ${outputdir}/${outputprefix}.$singlabel.UxxU.rm.bam.stat | wc -l )
singUxxUuniqnr=$(cat ${outputdir}/${outputprefix}.$singlabel.UxxU.rm.nr.bam.stat | wc -l | xargs printf "%'.f")
echo -e "single_UxxU_reads_dedup\t"${singUxxUuniqnr} >> ${out_file}
singUxxUuniqnr=$(cat ${outputdir}/${outputprefix}.$singlabel.UxxU.rm.nr.bam.stat | wc -l )

echo "" >> ${out_file}
nonelabeltotal=$(cat ${outputdir}/${outputprefix}.$nonelabel.sam.gz.stat | wc -l | xargs printf "%'.f")
echo -e "none_reads\t"${nonelabeltotal} >> ${out_file}
noneUUuniq=$(cat ${outputdir}/${outputprefix}.$nonelabel.UU.rm.bam.stat | wc -l | xargs printf "%'.f")
noneUxxUuniq=$(cat ${outputdir}/${outputprefix}.$nonelabel.UxxU.rm.bam.stat | wc -l | xargs printf "%'.f")
echo -e "none_UU_reads\t"${noneUUuniq} >> ${out_file}
echo -e "none_UxxU_reads\t"${noneUxxUuniq} >> ${out_file}
noneUUuniq=$(cat ${outputdir}/${outputprefix}.$nonelabel.UU.rm.bam.stat | wc -l )
noneUxxUuniq=$(cat ${outputdir}/${outputprefix}.$nonelabel.UxxU.rm.bam.stat | wc -l )
noneUUuniqnr=$(cat ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.bam.stat| wc -l | xargs printf "%'.f")
noneUxxUuniqnr=$(cat ${outputdir}/${outputprefix}.$nonelabel.UxxU.rm.nr.bam.stat | wc -l | xargs printf "%'.f")
echo -e "none_UU_reads_dedup\t"${noneUUuniqnr} >> ${out_file}
echo -e "none_UxxU_reads_dedup\t"${noneUxxUuniqnr} >> ${out_file}
noneUUuniqnr=$(cat ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.bam.stat| wc -l )
noneUxxUuniqnr=$(cat ${outputdir}/${outputprefix}.$nonelabel.UxxU.rm.nr.bam.stat | wc -l )
none_self_lig=$( cat ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.pet.span.xls | grep "second/best<0.95" -A5 |  awk -F '[\t]' '{if(NR==4)print $2}' )
none_self_lig=$( printf "%'.f" ${none_self_lig} )
echo -e "none_Self-ligation_PET\t"${none_self_lig} >> ${out_file}
none_intra_chr_pet=$( cat ${outputdir}/${outputprefix}.$nonelabel.UU.rm.nr.pet.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==5)print $2}' )
none_inter_chr_pet=$( cat ${outputdir}/${outputprefix}.$pairlabel.UU.rm.nr.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==2)print $2}' )
none_inter_lig_all=$( echo "${none_intra_chr_pet} + ${none_inter_chr_pet}" | bc )
none_inter_lig_all=$( printf "%'.f" ${none_inter_lig_all} )
echo -e "none_Inter-ligation_PET\t"${none_inter_lig_all} >> ${out_file}

none_intra_chr_pet=$( printf "%'.f" ${none_intra_chr_pet} )
echo -e "none_Intra-chr_PET\t"${none_intra_chr_pet} >> ${out_file}
none_inter_chr_pet=$( printf "%'.f" ${none_inter_chr_pet} )
echo -e "none_Inter-chr_PET\t"${none_inter_chr_pet} >> ${out_file}

echo "" >> ${out_file}
redun=$( echo "(${pairUUuniq}+${pairUxxUuniq}+${singUxxUuniq}+${noneUUuniq}+${noneUxxUuniq} - ${pairUUuniqnr}-${pairUxxUuniqnr}-${singUxxUuniqnr}-${noneUUuniqnr}-${noneUxxUuniqnr}) / (${pairUUuniq}+${pairUxxUuniq}+${singUxxUuniq}+${noneUUuniq}+${noneUxxUuniq})" | bc -l )
redun=$( printf %.2f ${redun} )
echo -e "Redundancy\t"${redun} >> ${out_file}


np=$( cat ${outputdir}/${outputprefix}_peaks.narrowPeak | wc -l )
n_peak=$( printf "%'.f" ${np} )
echo -e "Peak\t"$n_peak >> ${out_file}

self_lig=$( cat ${outputdir}/${outputprefix}.pet.span.xls | grep "second/best<0.95" -A5 |  awk -F '[\t]' '{if(NR==4)print $2}' )
self_lig=$( printf "%'.f" ${self_lig} )
echo -e "Self-ligation_PET\t"${self_lig} >> ${out_file}

intra_chr_pet=$( cat ${outputdir}/${outputprefix}.pet.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==5)print $2}' )
inter_chr_pet=$( cat ${outputdir}/${outputprefix}.pet.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==2)print $2}' )
pet_ratio=$( echo "${intra_chr_pet} / ${inter_chr_pet}" | bc -l )
pet_ratio=$( printf %.2f ${pet_ratio} )
inter_lig_all=$( echo "${intra_chr_pet} + ${inter_chr_pet}" | bc )
inter_lig_all=$( printf "%'.f" ${inter_lig_all} )
echo -e "Inter-ligation_PET\t"${inter_lig_all} >> ${out_file}

intra_chr_pet=$( printf "%'.f" ${intra_chr_pet} )
echo -e "Intra-chr_PET\t"${intra_chr_pet} >> ${out_file}
inter_chr_pet=$( printf "%'.f" ${inter_chr_pet} )
echo -e "Inter-chr_PET\t"${inter_chr_pet} >> ${out_file}


echo -e "ratio_of_intra/inter_PET\t"${pet_ratio} >> ${out_file}

echo "" >> ${out_file}
echo "for cluster" >> ${out_file}
total_cluster_number=$(zcat ${outputdir}/${outputprefix}*chiasig.gz | awk '$7 != 1{print}' | wc -l)
total_cluster_number=$( printf "%'.f" ${total_cluster_number} )
echo -e "PET_cluster\t"${total_cluster_number} >> ${out_file}

intra_cluster=$( zcat ${outputdir}/${outputprefix}*cis.chiasig.gz | awk '$7 >=2 {print}' | wc -l )
inter_cluster=$( zcat ${outputdir}/${outputprefix}*trans.chiasig.gz | awk '$7 >=2 {print}' | wc -l)
intra_cluster=$( printf "%'.f" ${intra_cluster} )
echo -e "Intra-chr_PET_cluster\t"${intra_cluster} >> ${out_file}
for i in $(seq 2 10)
do
    intra_pets_number=$(zcat ${outputdir}/${outputprefix}*cis.chiasig.gz | awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | xargs printf "%'.f")
    echo -e "pets_number_"${i}"\t"${intra_pets_number} >> ${out_file}
done
echo -e "pets_number>10\t"$(zcat *cis.chiasig.gz | awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

inter_cluster=$( printf "%'.f" ${inter_cluster} )
echo -e "Inter-chr_PET_cluster\t"${inter_cluster} >> ${out_file}
for i in $(seq 2 10)
do
    inter_pets_number=$(zcat ${outputdir}/${outputprefix}*trans.chiasig.gz | awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | xargs printf "%'.f")
    echo -e "pets_number_"${i}"\t"${inter_pets_number} >> ${out_file}
done
echo -e "pets_number>10\t"$(zcat *trans.chiasig.gz | awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

ns=$( grep -v : ${outputdir}/${outputprefix}.e${extbp}.clusters.cis.chiasig.BE2.sigf.interactions | wc -l )
ns=$( printf "%'.f" ${ns} )
echo -e "Significant_interaction>=${PET} \t"$ns >> ${out_file}

echo "" >> ${out_file}
echo "for cluster by given anchor" >> ${out_file}
total_cluster_number=$(cat ${outputdir}/${outputprefix}.forcluster.ipet.loops | awk '$7 != 1{print}' | wc -l)
total_cluster_number=$( printf "%'.f" ${total_cluster_number} )
echo -e "PET_cluster\t"${total_cluster_number} >> ${out_file}

intra_cluster=$( cat ${outputdir}/${outputprefix}.forcluster.ipet.loops | awk '$1==$4' | awk '$7 >=2 {print}' | wc -l )
inter_cluster=$( cat ${outputdir}/${outputprefix}.forcluster.ipet.loops | awk '$1!=$4' | awk '$7 >=2 {print}' | wc -l)
intra_cluster=$( printf "%'.f" ${intra_cluster} )
echo -e "Intra-chr_PET_cluster\t"${intra_cluster} >> ${out_file}
for i in $(seq 2 10)
do
    intra_pets_number=$( cat ${outputdir}/${outputprefix}.forcluster.ipet.loops | awk '$1==$4' | awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | xargs printf "%'.f")
    echo -e "pets_number_"${i}"\t"${intra_pets_number} >> ${out_file}
done
echo -e "pets_number>10\t"$( cat ${outputdir}/${outputprefix}.forcluster.ipet.loops | awk '$1==$4' | awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

inter_cluster=$( printf "%'.f" ${inter_cluster} )
echo -e "Inter-chr_PET_cluster\t"${inter_cluster} >> ${out_file}
for i in $(seq 2 10)
do
    inter_pets_number=$( cat ${outputdir}/${outputprefix}.forcluster.ipet.loops | awk '$1!=$4' | awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | xargs printf "%'.f")
    echo -e "pets_number_"${i}"\t"${inter_pets_number} >> ${out_file}
done
echo -e "pets_number>10\t"$( cat ${outputdir}/${outputprefix}.forcluster.ipet.loops | awk '$1!=$4' | awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

ns=$( cat ${outputdir}/${outputprefix}.forcluster.ipet.loops | awk '$7 != 1{print}' | awk '$10+0.0 > 0.001' | wc -l )
ns=$( printf "%'.f" ${ns} )
echo -e "Significant_interaction>=${PET} \t"$ns >> ${out_file}
