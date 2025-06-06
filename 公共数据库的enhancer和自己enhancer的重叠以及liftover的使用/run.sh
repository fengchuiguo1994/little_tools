grep ">" forebrain.txt | awk -F"|" '{print $2}' | awk -F":" -v OFS="\t" '{print $1,$2}' | awk -F"-" -v OFS="\t" '{print $1,$2}' > forebrain.bed
grep ">" hindbrain.txt | awk -F"|" '{print $2}' | awk -F":" -v OFS="\t" '{print $1,$2}' | awk -F"-" -v OFS="\t" '{print $1,$2}' > hindbrain.bed

CrossMap.py bed --unmap-file forebrain.umap.bed ~/Tools/CrossMap/mm9ToMm10.over.chain.gz forebrain.bed forebrain.mm10.bed
CrossMap.py bed --unmap-file hindbrain.umap.bed ~/Tools/CrossMap/mm9ToMm10.over.chain.gz hindbrain.bed hindbrain.mm10.bed


sed '1d' S.txt | awk -F":" -v OFS="\t" '{print $1,$2}' | awk -F"-" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2-10000,$3+10000,$1":"$2"-"$3}'> mousebrain.bed


awk -v OFS="\t" '{print $1,$2,$3}' hindbrain.bed | intersectBed -a mousebrain.bed -b - -wao > mousebrain.hindbrain.bed
awk -v OFS="\t" '{print $1,$2,$3}' forebrain.bed | intersectBed -a mousebrain.bed -b - -wao > mousebrain.forebrain.bed

