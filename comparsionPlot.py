# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 14:50:04 2023

@author: laniq
"""
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from scipy.stats import linregress

def displayComparisionPlots(displacementArray_x, displacementArray_y, movingMaskFileName, save):
    """
    displayComparisionPlots displays two different sets of displacement data to one another 

    :param displacementArray_x, displacementArray_y: first set of displacement data in the x and y fields
    :param movingMaskFileName: file name of the mask to apply to image
    :param save: whether or not to save the resulting displayed image
    """ 
    pix_size = 0.908
    beat = 1
    
    #Emma's data
    pos = {}
    pos['x'] = np.loadtxt('fibrotug_2/comparisionPlots' + '/beat%i' % beat + '_col.txt')
    pos['y'] = np.loadtxt('fibrotug_2/comparisionPlots' + '/beat%i' % beat + '_row.txt')
    
    #Elastix data
    elastix_x = displacementArray_x.T
    elastix_y = (displacementArray_y).T
    
    mask = np.asarray(imread(movingMaskFileName, as_gray=True)).T
    
    #use these values to find coordinates to track on elastix's displacement field
    x_0 = np.round(pos['x'][:,0],0).astype(int)
    y_0 = np.round(pos['y'][:,0],0).astype(int)
    x_25 = np.round(pos['x'][:,25],0).astype(int)
    y_25 = np.round(pos['y'][:,25],0).astype(int)
    
    x_0_float = pos['x'][:,0]
    y_0_float = pos['y'][:,0]
    x_25_float = pos['x'][:,25]
    y_25_float = pos['y'][:,25]
    
    #remove any coordinates that are not in elastix's displacement fields
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
    
    #Calculate the displacement based on original and contracted images
    disp = {}
    disp['x'] = x_25_float - x_0_float
    disp['y'] = y_25_float - y_0_float
    
    elastix_disp = {}
    elastix_disp['x'] = elastix_x[x_25, y_25]
    elastix_disp['y'] = elastix_y[x_25, y_25]
    
    #linear regression
    lres_x = linregress(disp['x'],  elastix_disp['x'])
    slope_x = lres_x.slope
    intercept_x = lres_x.intercept
    
    lres_y = linregress(disp['y'],  elastix_disp['y'])
    slope_y = lres_y.slope
    intercept_y = lres_y.intercept
    
    #r^2
    r2_x = linregress(disp['x'],  elastix_disp['x']).rvalue**2
    r2_y = linregress(disp['y'],  elastix_disp['y']).rvalue**2
    
    #Display images
    fig, axs = plt.subplots(1, 2)
    axs[0].scatter(disp['x'], elastix_disp['x'], marker = '.', color='gold')
    x = [-12,12]
    axs[0].plot(x, x, color = 'blue')
    axs[0].set_title('x direction')
    axs[0].annotate('R2=%2.3f' % r2_x, (0.05,0.9), xycoords='axes fraction')
    axs[0].annotate('y=%2.3f' % slope_x + 'x + %2.3f' % intercept_x, (0.05,0.85), xycoords='axes fraction')
    axs[0].set_xlabel("emmas")
    axs[0].set_ylabel("elastix")
    
    axs[1].scatter(disp['y'], elastix_disp['y'], marker = '.', color='gold')
    x2 = [-3,5]
    axs[1].plot(x2, x2, color = 'blue')
    axs[1].set_title('y direction')
    axs[1].annotate('R2=%2.3f' % r2_y, (0.05,0.9), xycoords='axes fraction')
    axs[1].annotate('y=%2.3f' % slope_y + 'x + %2.3f' % intercept_y, (0.05,0.85), xycoords='axes fraction')
    axs[1].set_xlabel("emmas")
    
    fig.suptitle('FibroTug Elastix v Emmas Comparsion Plots', fontsize=12)
    
    if save:
        plt.savefig('fibrotug_2/comparisionPlots/comparison_plot.png', dpi = 180)


def scatter(displacementArray_y, movingMaskFileName):
    beat = 1
    
    #Emma's data
    pos = {}
    pos['x'] = np.loadtxt('fibrotug_2/comparisionPlots' + '/beat%i' % beat + '_col.txt')
    pos['y'] = np.loadtxt('fibrotug_2/comparisionPlots' + '/beat%i' % beat + '_row.txt')
    
    mask = np.asarray(imread(movingMaskFileName, as_gray=True)).T
    
    #use these values to find coordinates to track on elastix's displacement field
    x_0 = np.round(pos['x'][:,0],0).astype(int)
    y_0 = np.round(pos['y'][:,0],0).astype(int)
    x_25 = np.round(pos['x'][:,25],0).astype(int)
    y_25 = np.round(pos['y'][:,25],0).astype(int)
    
    x_0_float = pos['x'][:,0]
    y_0_float = pos['y'][:,0]
    x_25_float = pos['x'][:,25]
    y_25_float = pos['y'][:,25]
    
    #remove any coordinates that are not in elastix's displacement fields
    # x_0_float[mask[x_25, y_25] == 0] = -1
    # y_0_float[mask[x_25, y_25] == 0] = -1
    # x_25_float[mask[x_25, y_25] == 0] = -1
    # y_25_float[mask[x_25, y_25] == 0] = -1
    
    # x_0_float[mask[x_25, y_25] == 0] = -1
    # y_0_float[mask[x_25, y_25] == 0] = -1
    # x_25_float[mask[x_25, y_25] == 0] = -1
    # y_25_float[mask[x_25, y_25] == 0] = -1
    
    # x_0_float = x_0_float[x_0_float >-1]
    # y_0_float = y_0_float[y_0_float > -1]
    # x_25_float = x_25_float[x_25_float > -1]
    # y_25_float = y_25_float[y_25_float > -1]
    
    
    plt.figure()
    im = plt.imshow(displacementArray_y, vmin=-5, vmax=5)
    plt.scatter(x_0_float, y_0_float, c=y_25_float - y_0_float, vmin=-5, vmax=5, edgecolors='k', linewidth=0.5)
    plt.colorbar(im)
    plt.axis('off')
    plt.title('Displacement Y')





