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
beat = 0

pos = {}
pos['x'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_col.txt')
pos['y'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_row.txt')

analytical_x = np.load('comparisionPlots/analytical_dispx_45.npy').T
analytical_y = (np.load('comparisionPlots/analytical_dispy_45.npy')*-1).T

mask = np.asarray(imread('transformingTestImage/Mask_41.png', as_gray=True)).T

x_0 = np.round(pos['x'][:,0],0).astype(int)
y_0 = np.round(pos['y'][:,0],0).astype(int)
x_41= np.round(pos['x'][:,41],0).astype(int)
y_41 = np.round(pos['y'][:,41],0).astype(int)

x_0_float = pos['x'][:,0]
y_0_float = pos['y'][:,0]
x_41_float = pos['x'][:,41]
y_41_float = pos['y'][:,41]

x_0_float[mask[x_41, y_41] == 0] = -1
y_0_float[mask[x_41, y_41] == 0] = -1
x_41_float[mask[x_41, y_41] == 0] = -1
y_41_float[mask[x_41, y_41] == 0] = -1

x_0_float = x_0_float[x_0_float >-1]
y_0_float = y_0_float[y_0_float > -1]
x_41_float = x_41_float[x_41_float > -1]
y_41_float = y_41_float[y_41_float > -1]

x_0[mask[x_41, y_41] == 0] = -1 
y_0[mask[x_41, y_41] == 0] = -1 
x_41[mask[x_41, y_41] == 0] = -1 
y_41[mask[x_41, y_41] == 0] = -1 
 
x_0 = x_0[x_0 > -1] 
y_0 = y_0[y_0 > -1] 
x_41 = x_41[x_41 > -1] 
y_41 = y_41[y_41 > -1] 

disp = {}
disp['x'] = x_41_float - x_0_float
disp['y'] = y_41_float - y_0_float

analytical_disp = {}
analytical_disp['x'] = analytical_x[x_41, y_41]
analytical_disp['y'] = analytical_y[x_41, y_41]

r2_x = linregress(disp['x'],  analytical_disp['x']).rvalue**2
r2_y = linregress(disp['y'],  analytical_disp['y']).rvalue**2

fig, axs = plt.subplots(1, 2)
axs[0].scatter(analytical_disp['x'], disp['x'], marker = '.', color='gold')
x = [-10,10]
axs[0].plot(x, x, color = 'blue')
axs[0].set_title('x direction')
axs[0].annotate('R2=%2.3f' % r2_x, (0.05,0.9), xycoords='axes fraction')
axs[0].set_ylabel("emmas")
axs[0].set_xlabel("analytical")

axs[1].scatter(analytical_disp['y'], disp['y'], marker = '.', color='gold')
x2 = [-8,2]
axs[1].plot(x2, x2, color = 'blue')
axs[1].set_title('y direction')
axs[1].annotate('R2=%2.3f' % r2_y, (0.05,0.9), xycoords='axes fraction')
axs[1].set_xlabel("analytical")

plt.savefig('comparisionPlots/Emma_v_Analytical_plot.png')







