# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 20:00:24 2023

@author: laniq
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from scipy.stats import linregress
from mpl_toolkits.axes_grid1 import make_axes_locatable
pix_size = 0.908
beat = 0

pos = {}
pos['x'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_col.txt')
pos['y'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_row.txt')

analytical_x = np.load('comparisionPlots/analytical_dispx_45.npy')
analytical_y = (np.load('comparisionPlots/analytical_dispy_45.npy')*-1)

elastix_x = np.load('comparisionPlots/displacement_x.npy').T
elastix_y = (np.load('comparisionPlots/displacement_y.npy')).T

mask = np.asarray(imread('transformingTestImage/Mask_41.png', as_gray=True))

analytical_disp = {}
analytical_disp['x'] = analytical_x[:,:]
analytical_disp['y'] = analytical_y[:,:]

analytical_disp['y'][mask==0] = np.nan
analytical_disp['x'][mask==0] = np.nan


#Plot images
fig, axs = plt.subplots(1, 2, sharey=True, figsize=[30,30])
im3 = axs[1].imshow(analytical_disp['x'], vmin = -10, vmax = 10, cmap='BrBG')

im2 = axs[0].imshow(analytical_disp['y'], vmin = -3, vmax = 7, cmap = 'RdYlBu')

axs[0].axis('off')
axs[1].axis('off')
axs[0].set_title('Displacement Field Y', fontsize=30)
axs[1].set_title('Displacement Field X', fontsize=30)

plt.savefig('analytical_disp_insilico', dpi = 200)