# python
### 安装
```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
```

### conda
安装
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
环境配置
```
conda create -n bioTools python=3.8 # 安装某个环境
conda info -e # 列出已经安装的环境
source activate bioTools # 切入某个环境
pip install deeptools
pip install macs2
pip install cutadapt
conda install -c bioconda perl-app-cpanminus
conda install -c r r-base=3.6
conda install -c bioconda bwa
conda install -c bioconda samtools=1.9
conda install -c bioconda fastp
```