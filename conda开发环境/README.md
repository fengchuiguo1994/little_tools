```
conda env export > 10.73.29.10.45777.environment.yml
conda env export -n CellDART > CellDART_environment.yml


conda env list
for i in GDAL GDALS4 R4 R4bak RNAfold Rarchr SpatialGlue2 alfred cell2loc_env deeptools fastp  hicpeaks kmc ldsc meme neuronVis py38 py3hapcut2 snakemake st starfusion; do conda env export -n ${i} > ${i}_environment.yml; done


```