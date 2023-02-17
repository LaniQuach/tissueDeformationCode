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

elastix_x = np.load('comparisionPlots/displacement_x.npy').T
elastix_y = (np.load('comparisionPlots/displacement_y.npy')).T

mask = np.asarray(imread('transformingTestImage/Mask_41.png', as_gray=True)).T

x_0 = np.round(pos['x'][:,0],0).astype(int)
y_0 = np.round(pos['y'][:,0],0).astype(int)
x_41 = np.round(pos['x'][:,41],0).astype(int)
y_41 = np.round(pos['y'][:,41],0).astype(int)

x_0[mask[x_41, y_41] == 0] = -1 
y_0[mask[x_41, y_41] == 0] = -1 
x_41[mask[x_41, y_41] == 0] = -1 
y_41[mask[x_41, y_41] == 0] = -1 
 
x_0 = x_0[x_0 > -1] 
y_0 = y_0[y_0 > -1] 
x_41 = x_41[x_41 > -1] 
y_41 = y_41[y_41 > -1] 

elastix_disp = {}
elastix_disp['x'] = elastix_x[x_41, y_41]
elastix_disp['y'] = elastix_y[x_41, y_41]

analytical_disp = {}
analytical_disp['x'] = analytical_x[x_41, y_41]
analytical_disp['y'] = analytical_y[x_41, y_41]

#r^2 values
r2_x = linregress(analytical_disp['x'],  elastix_disp['x']).rvalue**2
r2_y = linregress(analytical_disp['y'],  elastix_disp['y']).rvalue**2

#error
disp_analytical = np.vstack([analytical_disp['x'], analytical_disp['y']]).T
disp_elastix = np.vstack([elastix_disp['x'], elastix_disp['y']]).T
error_elastix = np.linalg.norm(disp_analytical-disp_elastix)

#linear regression
lres_x = linregress(analytical_disp['x'],  elastix_disp['x'])
slope_x = lres_x.slope
intercept_x = lres_x.intercept

lres_y = linregress(analytical_disp['y'],  elastix_disp['y'])
slope_y = lres_y.slope
intercept_y = lres_y.intercept

#plotting
fig, axs = plt.subplots(1, 2)
axs[0].scatter(analytical_disp['x'], elastix_disp['x'], marker = '.', color='gold')
x = [-10,10]
axs[0].plot(x, x, color = 'blue')
axs[0].set_title('x direction')
axs[0].annotate('R2=%2.3f' % r2_x, (0.05,0.9), xycoords='axes fraction')
axs[0].annotate('y=%2.3f' % slope_x + 'x + %2.3f' % intercept_x, (0.05,0.85), xycoords='axes fraction')
axs[0].set_xlabel("analytical")
axs[0].set_ylabel("elastix")

axs[1].scatter(analytical_disp['y'], elastix_disp['y'], marker = '.', color='gold')
x2 = [-8,2]
axs[1].plot(x2, x2, color = 'blue')
axs[1].set_title('y direction')
axs[1].annotate('R2=%2.3f' % r2_y, (0.05,0.9), xycoords='axes fraction')
axs[1].annotate('y=%2.3f' % slope_y + 'x + %2.3f' % intercept_y, (0.05,0.85), xycoords='axes fraction')
axs[1].set_xlabel("analytical")

fig.suptitle(' In-Silico Elastix v Analytical Comparsion Plots', fontsize=12)
txt="error: %2.3f" % error_elastix
plt.figtext(0.5, -0.01, txt, wrap=True, horizontalalignment='center', fontsize=9)

plt.savefig('comparisionPlots/Analytical_v_Elastics_plot.png')







