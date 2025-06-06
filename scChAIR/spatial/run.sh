# sbatch -p ruan_cpu -J C02533C1.run -o C02533C1.run.out -e C02533C1.run.err -N 1 -n 30 --mem=550G --wrap=" bash C02533C1.sh "

# sbatch -p ruan_cpu -J A02677B5.run -o A02677B5.run.out -e A02677B5.run.err -N 1 -n 30 --mem=550G --wrap=" bash A02677B5.sh "

# singularity exec /data/home/ruanlab/huangxingyu/Tools/SAW_v7.1.sif spatialCluster -i /data/home/ruanlab/huangxingyu/Haoxi20230215/spatial/DemoData/C02533C1result/result/04.tissuecut/C02533C1.tissue.gef -s 50 -r 1.0 -o /data/home/ruanlab/huangxingyu/Haoxi20230215/spatial/DemoData/C02533C1result/result/05.spatialcluster/C02533C1.bin50_1.0.spatial.cluster.h5ad
