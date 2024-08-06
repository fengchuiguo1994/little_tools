awk '{len=int(($3-$2)/10);print $1"\t"$2-len"\t"$3+len}' ~/multi-ChIA/meizheng/Drosophila/ChIPSeq/RNAPII/dm3.narrow_peaks.narrowPeak | grep -v track > dm3.RNAPII.anchor

java -Xmx10G -cp /public/home/xyhuang/longread_pipeline/longread_pipeline/program/LGL.jar LGL.chiapet.PetClusterWithGivenAnchors ../dm3.RNAPII.pre_cluster.sorted dm3.RNAPII.anchor dm3.RNAPII.cluster 1
awk '{if($13>=2)print}' < dm3.RNAPII.cluster.cluster.txt > dm3.RNAPII.cluster.filtered
cut -f1-3,7-9,13-15 < dm3.RNAPII.cluster.filtered > dm3.RNAPII.cluster.txt

###### calculation of p-values
cut -f1-3 < ../dm3.RNAPII.ipet > dm3.RNAPII.aln
cut -f4-6 < ../dm3.RNAPII.ipet >> dm3.RNAPII.aln
cut -f1-3,13 < dm3.RNAPII.cluster.filtered > dm3.RNAPII.cluster.filtered.anchor1
cut -f7-9,13 < dm3.RNAPII.cluster.filtered > dm3.RNAPII.cluster.filtered.anchor2

for y in anchor1 anchor2
do
    java -Xmx10G -cp /public/home/xyhuang/longread_pipeline/longread_pipeline/program/LGL.jar LGL.shortReads.TagCountInGivenRegions dm3.RNAPII.aln dm3.RNAPII.cluster.filtered.${y} dm3.RNAPII.cluster.filtered.${y}.tagCount.txt 1 2
done

## generate the global tag count for global density
wc -l dm3.RNAPII.aln | sed 's/ /\t/g' |  cut -f1 > dm3.RNAPII.nTags.txt
## calculate p-value
cp dm3.RNAPII.nTags.txt nTags.txt
paste dm3.RNAPII.cluster.filtered.anchor1.tagCount.txt  dm3.RNAPII.cluster.filtered.anchor2.tagCount.txt |  cut -f4,5,10  > dm3.RNAPII.petCount.tagCount.txt
cp dm3.RNAPII.petCount.tagCount.txt data.txt
R --vanilla < /public/home/xyhuang/longread_pipeline/longread_pipeline/program/hypergeometric.r

mv -f result.txt dm3.RNAPII.pvalue.hypergeo.txt
paste dm3.RNAPII.cluster.filtered dm3.RNAPII.petCount.tagCount.txt dm3.RNAPII.pvalue.hypergeo.txt | cut -f1-3,7-9,13-15,19-24 > dm3.RNAPII.cluster.withpvalue.txt
awk -v pvalue_cutoff=0.05 '{if($13<pvalue_cutoff+0.0)print}' < dm3.RNAPII.cluster.withpvalue.txt > dm3.RNAPII.cluster.FDRfiltered.txt

awk '{print $7"\t"$8}' dm3.RNAPII.cluster.FDRfiltered.txt | LANG=C sort | uniq -c | awk '{if($2>=10){a+=$1;b+=$1*$3;f+=$1;g+=$1*$3}else{c[$2]+=$1;d[$2]+=$1*$3;f+=$1;g+=$1*$3}}END{for(i in c){print c[i]"\t"i"\t"d[i]};print a"\t10\t"b;print f"\t11\t"g}' | LANG=C sort -k2,2n | awk 'BEGIN{print "PET counts\tNo. of clusters\tNo.intra-chrom clusters\tNo.inter-chrom clusters\tPercent of intra-chrom clusters"}{if($2==10){print ">="$2"\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}else if($2==11){print "Total\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}else{print $2"\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}}' > dm3.RNAPII.PET_count_distribution.txt


awk '{if ($8==1){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$12"\t"$13}}' dm3.RNAPII.cluster.FDRfiltered.txt > dm3.RNAPII.giveanchor.cluster.curve

#awk '{if($3=="-"){print $1"\t"$2"\t"$5}else{print $1"\t"$5"\t"$2}}' ../dm3.RNAPII.spet | sort -k1,1 -k2,2n > dm3.RNAPII.coverage
#cat dm3.RNAPII.coverage | ~/tools/bedItemOverlapCount -chromSize=/public/home/rqzheng/tools/longread_pipeline/program/rapeseed.chromSize.txt test stdin > B2CP-007.bedGraph
#~/tools/bedGraphToBigWig B2CP-007.bedGraph /public/home/rqzheng/tools/longread_pipeline/program/rapeseed.chromSize.txt B2CP-007.RNAP2.coverage.bigWig


awk '$7>4' dm3.RNAPII.cluster.FDRfiltered.txt > dm3.RNAPII.cluster.FDRfiltered.gt4.txt
python ~/longread_pipeline/longread_pipeline/program/split_region.py dm3.RNAPII.cluster.FDRfiltered.gt4.txt dm3.RNAPII.cluster.phase.gt4
