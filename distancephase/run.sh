for line in {2..9}
do
awk -v pet=$line '{if($7==pet && $8==1){print $9}}' ../A549.CTCF.cluster.FDRfiltered.txt > pet$line.txt
done
awk -v pet=10 '{if($7>=pet && $8==1){print $9}}' ../A549.CTCF.cluster.FDRfiltered.txt > pet10.txt

awk -v pet=10 '{if($7>=pet){$7=pet;}if($8==1){print $7"\t"$9}}' ../A549.CTCF.cluster.FDRfiltered.txt > pet.txt

Rscript plot.r
# pick  petcount > 4


awk '$7>4' ../A549.CTCF.cluster.FDRfiltered.txt > A549.CTCF.cluster.FDRfiltered.gt4.txt

