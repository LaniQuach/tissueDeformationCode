# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 14:50:04 2023

@author: laniq
"""
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from scipy.stats import linregress

pix_size = 0.908
beat = 1

pos = {}
pos['x'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_col.txt')
pos['y'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_row.txt')

elastix_x = np.load('comparisionPlots/displacement_x.npy').T
elastix_y = (np.load('comparisionPlots/displacement_y.npy')*-1).T

mask = np.asarray(imread('transformingTestImage/Mask_25.png', as_gray=True)).T

x_0 = np.round(pos['x'][:,0],0).astype(int)
y_0 = np.round(pos['y'][:,0],0).astype(int)
x_25 = np.round(pos['x'][:,25],0).astype(int)
y_25 = np.round(pos['y'][:,25],0).astype(int)

x_0_float = pos['x'][:,0]
y_0_float = pos['y'][:,0]
x_25_float = pos['x'][:,25]
y_25_float = pos['y'][:,25]

x_0_float[mask[x_25, y_25] == 0] = -1
y_0_float[mask[x_25, y_25] == 0] = -1
x_25_float[mask[x_25, y_25] == 0] = -1
y_25_float[mask[x_25, y_25] == 0] = -1

x_0_float[mask[x_25, y_25] == 0] = -1
y_0_float[mask[x_25, y_25] == 0] = -1
x_25_float[mask[x_25, y_25] == 0] = -1
y_25_float[mask[x_25, y_25] == 0] = -1

x_0_float = x_0_float[x_0_float >-1]
y_0_float = y_0_float[y_0_float > -1]
x_25_float = x_25_float[x_25_float > -1]
y_25_float = y_25_float[y_25_float > -1]

x_0[mask[x_25, y_25] == 0] = -1 
y_0[mask[x_25, y_25] == 0] = -1 
x_25[mask[x_25, y_25] == 0] = -1 
y_25[mask[x_25, y_25] == 0] = -1 
 
x_0 = x_0[x_0 > -1] 
y_0 = y_0[y_0 > -1] 
x_25 = x_25[x_25 > -1] 
y_25 = y_25[y_25 > -1] 

disp = {}
disp['x'] = x_25_float - x_0_float
disp['y'] = y_25_float - y_0_float

elastix_disp = {}
elastix_disp['x'] = elastix_x[x_25, y_25]
elastix_disp['y'] = elastix_y[x_25, y_25]

r2_x = linregress(disp['x'],  elastix_disp['x']).rvalue**2
r2_y = linregress(disp['y'],  elastix_disp['y']).rvalue**2

fig, axs = plt.subplots(1, 2)
axs[0].scatter(disp['x'], elastix_disp['x'], marker = '.', color='gold')
x = [-12,12]
axs[0].plot(x, x, color = 'blue')
axs[0].set_title('x direction')
axs[0].annotate('R2=%2.3f' % r2_x, (0.05,0.9), xycoords='axes fraction')
axs[0].set_xlabel("emmas")
axs[0].set_ylabel("elastix")

axs[1].scatter(disp['y'], elastix_disp['y'], marker = '.', color='gold')
x2 = [-3,5]
axs[1].plot(x2, x2, color = 'blue')
axs[1].set_title('y direction')
axs[1].annotate('R2=%2.3f' % r2_y, (0.05,0.9), xycoords='axes fraction')
axs[1].set_xlabel("emmas")

plt.savefig('comparisionPlots/Emma_v_Elastics_plot.png')







