# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 08:34:29 2022

@author: Lani
"""

from skimage import io
import os
import matplotlib.pyplot as plt

# Load .tif stack
# path = 'C:\Users\laniq\OneDrive\Documents\cardio\tissueDeformationCode\fibertug\data\fibers_2'
tissue = 'A_0.41_01-14-2_rc_ds'
im = io.imread(tissue + '.tif')

# To save a frame
# path = 'C:\Users\laniq\OneDrive\Documents\cardio\tissueDeformationCode\fibertug\data\fibers_2'
frame = 100
ext = '.png'
# io.imsave(completeName, im[frame])

# To save all frames (it will create a bunch of files, so might want to create a specific folder for this)
ext = '.png'
for frame in range(0, im.shape[0]):
    completeName = tissue + '_%i' % frame + ext
    # completeName = os.path.join(path, filename)
    io.imsave(completeName, im[frame])
