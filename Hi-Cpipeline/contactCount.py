import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
from multiprocessing import Pool
# import open2c libraries
import bioframe
import cooler
import cooltools

resolution = 1000
clr  = cooler.Cooler('normal.hic2cool.mcool::resolutions/1000')
cvd = cooltools.expected_cis(clr=clr, smooth=False, aggregate_smoothed=False, nproc=5)

f, ax = plt.subplots(1,1)
ax.loglog(cvd['dist'], cvd['count.sum'])
ax.set(xlabel='separation, bp', ylabel='IC contact frequency')
ax.set_aspect(1.0)
ax.grid(lw=0.5)
f.savefig("test.test.pdf")