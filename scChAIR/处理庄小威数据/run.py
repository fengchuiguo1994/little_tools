import scanpy as sc
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sc.settings.verbosity = 3
sc.logging.print_header()
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=120, facecolor='white')

adatas = sc.read_h5ad("c92c5b85-61d4-4228-b8fc-24ce4a757b59.h5ad")
adatas.var['highly_variable'].value_counts()
adatas.X # mydat@assays$RNA@data[1:5,1:5]
tt = adatas.raw # tt = sc.AnnData(adatas.raw.X, obs=adatas.obs, var=adatas.var)
print(tt.X) # mydat@assays$RNA@counts[1:10,1:5]
adatas.obs.to_csv('mouseCortex.metadata.txt', index=False, sep='\t') # 保存matedata



mydat = sc.AnnData(adatas.raw.X, obs=adatas.obs, var=adatas.var)
pdf = PdfPages('mouseCortex.exp.pdf')
plt.figure()
sc.pl.highest_expr_genes(mydat, n_top=20)
plt.tight_layout()
pdf.savefig()
plt.close()
pdf.close()

# sc.pp.filter_cells(adatas, min_genes=200)
# sc.pp.filter_genes(adatas, min_cells=3)