# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K27ac.bowtie2.bam ../S1.bowtie2.bam ../S2.bowtie2.bam && samtools index SiHa.H3K27ac.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K4me1.bowtie2.bam ../S3.bowtie2.bam ../S4.bowtie2.bam  && samtools index SiHa.H3K4me1.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.H3K4me3.bowtie2.bam ../S5.bowtie2.bam ../S6.bowtie2.bam && samtools index SiHa.H3K4me3.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SiHa.Sp1.bowtie2.bam ../S7.bowtie2.bam ../S8.bowtie2.bam && samtools index SiHa.Sp1.bowtie2.bam  "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K27ac.bowtie2.bam ../SP9.bowtie2.bam ../SP10.bowtie2.bam && samtools index SP.H3K27ac.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K4me1.bowtie2.bam ../SP11.bowtie2.bam ../SP12.bowtie2.bam && samtools index SP.H3K4me1.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.H3K4me3.bowtie2.bam ../SP13.bowtie2.bam ../SP14.bowtie2.bam && samtools index SP.H3K4me3.bowtie2.bam "
# bsub -q q2680v2 -J merge -n 5 -o merge.out -e merge.err -R span[hosts=1] " samtools merge -@ 5 SP.Sp1.bowtie2.bam ../SP15.bowtie2.bam ../SP16.bowtie2.bam && samtools index SP.Sp1.bowtie2.bam "

# cp ../Input1.bowtie2.bam .
# cp ../Input2.bowtie2.bam .
# cp ../Input3.bowtie2.bam .
# cp ../Input4.bowtie2.bam .
# samtools index Input1.bowtie2.bam 
# samtools index Input2.bowtie2.bam
# samtools index Input3.bowtie2.bam 
# samtools index Input4.bowtie2.bam

# for i in SiHa.H3K27ac SiHa.H3K4me1 SiHa.H3K4me3 SP.H3K27ac SP.H3K4me1 SP.H3K4me3 Input1 Input3
# do
#     echo $i
#     bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bowtie2.bam -o ${i}.bowtie2.split -p dm3_ "
# done
# for i in SiHa.Sp1 SP.Sp1 Input2 Input4
# do
#     echo $i
#     bsub -q q2680v2 -J split${i} -n 1 -o split${i}.out -e split${i}.err -R span[hosts=1] " split_bam.py -i ${i}.bowtie2.bam -o ${i}.bowtie2.split -p mm10_ "
# done


# cat *.*.bowtie2.split.report.txt Input?.bowtie2.split.report.txt | grep -v n_exogenous | awk '{print $0"\t"1000000/$NF}'
# SiHa.H3K27ac.bowtie2.bam	2335182	0	0	0	5487531	6355	71981379	5000087	0.199997
# SiHa.H3K4me1.bowtie2.bam	1962676	0	0	0	5359482	10684	60774539	9212725	0.108546
# SiHa.H3K4me3.bowtie2.bam	2419709	0	0	0	5716841	18138	67630941	13737989	0.0727909
# SiHa.Sp1.bowtie2.bam	2038129	0	0	0	14639663	35247	19492293	42508434	0.0235247
# SP.H3K27ac.bowtie2.bam	4677058	0	0	0	7737138	2630	76900176	2026070	0.493566
# SP.H3K4me1.bowtie2.bam	2987562	0	0	0	7146645	3472	67217391	3742868	0.267175
# SP.H3K4me3.bowtie2.bam	3437354	0	0	0	7698955	5612	78977338	5726787	0.174618
# SP.Sp1.bowtie2.bam	3212986	0	0	0	11795122	20767	55865179	15381730	0.0650122
# Input1.bowtie2.bam	4646411	0	0	0	1884924	2067	12376493	2358857	0.423934
# Input2.bowtie2.bam	1257990	0	0	0	3376013	6922	2322333	8553218	0.116915
# Input3.bowtie2.bam	3483991	0	0	0	1934904	1057	14033987	623397	1.60411
# Input4.bowtie2.bam	1583138	0	0	0	3289502	7975	4726775	7460630	0.134037



# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K27ac.bowtie2.bam -c Input1.bowtie2.bam --spikeIn --csf 0.423934 --tsf 0.199997 -o SiHa.H3K27ac.bowtie2.narrowpeak --bw --frip "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K4me1.bowtie2.bam -c Input1.bowtie2.bam --spikeIn --csf 0.423934 --tsf 0.108546 -o SiHa.H3K4me1.bowtie2.narrowpeak --bw --frip "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K4me3.bowtie2.bam -c Input1.bowtie2.bam --spikeIn --csf 0.423934 --tsf 0.0727909 -o SiHa.H3K4me3.bowtie2.narrowpeak --bw --frip "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.Sp1.bowtie2.bam -c Input2.bowtie2.bam --spikeIn --csf 0.116915 --tsf 0.0235247 -o SiHa.Sp1.bowtie2.narrowpeak --bw --frip "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K27ac.bowtie2.bam -c Input3.bowtie2.bam --spikeIn --csf 1.60411 --tsf 0.493566 -o SP.H3K27ac.bowtie2.narrowpeak --bw --frip "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K4me1.bowtie2.bam -c Input3.bowtie2.bam --spikeIn --csf 1.60411 --tsf 0.267175 -o SP.H3K4me1.bowtie2.narrowpeak --bw --frip "
# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K4me3.bowtie2.bam -c Input3.bowtie2.bam --spikeIn --csf 1.60411 --tsf 0.174618 -o SP.H3K4me3.bowtie2.narrowpeak --bw --frip "

# bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.Sp1.bowtie2.bam -c Input4.bowtie2.bam --spikeIn --csf 0.134037 --tsf 0.0650122 -o SP.Sp1.bowtie2.narrowpeak --bw --frip "




bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K27ac.bowtie2.bam -c Input1.bowtie2.bam --spikeIn --csf 0.423934 --tsf 0.199997 -o SiHa.H3K27ac.bowtie2.broadpeak --bw --frip --broad "
bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K4me1.bowtie2.bam -c Input1.bowtie2.bam --spikeIn --csf 0.423934 --tsf 0.108546 -o SiHa.H3K4me1.bowtie2.broadpeak --bw --frip --broad "
bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.H3K4me3.bowtie2.bam -c Input1.bowtie2.bam --spikeIn --csf 0.423934 --tsf 0.0727909 -o SiHa.H3K4me3.bowtie2.broadpeak --bw --frip --broad "

bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SiHa.Sp1.bowtie2.bam -c Input2.bowtie2.bam --spikeIn --csf 0.116915 --tsf 0.0235247 -o SiHa.Sp1.bowtie2.broadpeak --bw --frip --broad "

bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K27ac.bowtie2.bam -c Input3.bowtie2.bam --spikeIn --csf 1.60411 --tsf 0.493566 -o SP.H3K27ac.bowtie2.broadpeak --bw --frip --broad "
bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K4me1.bowtie2.bam -c Input3.bowtie2.bam --spikeIn --csf 1.60411 --tsf 0.267175 -o SP.H3K4me1.bowtie2.broadpeak --bw --frip --broad "
bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.H3K4me3.bowtie2.bam -c Input3.bowtie2.bam --spikeIn --csf 1.60411 --tsf 0.174618 -o SP.H3K4me3.bowtie2.broadpeak --bw --frip --broad "

bsub -q q2680v2 -J test -n 8 -o test.out -e test.err -R span[hosts=1] " spiker.py -t SP.Sp1.bowtie2.bam -c Input4.bowtie2.bam --spikeIn --csf 0.134037 --tsf 0.0650122 -o SP.Sp1.bowtie2.broadpeak --bw --frip --broad "
