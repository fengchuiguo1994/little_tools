# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:49:34 2019

@author: huang
"""
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

df = ["chr1",2902,4067,'chr1',24099,25515]
before = int((df[1]+df[2])/2)
after = int((df[4]+df[5])/2)
min2 = int((before+after)/2)
Path = mpath.Path
axs = plt.subplot(211)
pp1 = mpatches.PathPatch(
Path([(before, 0), (min2,2*2), (after, 0),(before, 0)],
     [Path.MOVETO, Path.CURVE3, Path.CURVE3,Path.CLOSEPOLY]),
fc="none", transform=axs.transData,color="red")
axs.add_patch(pp1)
axs.axis([0,35000,0,8])