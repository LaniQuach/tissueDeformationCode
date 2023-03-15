# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 14:50:04 2023

@author: laniq
"""
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from scipy.stats import linregress


def displayComparisionPlots(displacementArray_x, displacementArray_y, analyticalDisp_x, analyticalDisp_y, movingMaskFileName, save):
    """
    displayComparisionPlots displays two different sets of displacement data to one another 

    :param displacementArray_x, displacementArray_y: first set of displacement data in the x and y fields
    :param movingMaskFileName: file name of the mask to apply to image
    :param save: whether or not to save the resulting displayed image
    """ 
    pix_size = 0.908
    beat = 0
    
    #Emma's data
    pos = {}
    pos['x'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_col.txt')
    pos['y'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_row.txt')
    
    #Elastix data
    elastix_x = displacementArray_x.T
    elastix_y = (displacementArray_y).T
    
    analytical_x = analyticalDisp_x.T
    analytical_y = analyticalDisp_y.T
    
    mask = np.asarray(imread(movingMaskFileName, as_gray=True)).T
    
    #use these values to find coordinates to track on elastix's displacement field
    x_0 = np.round(pos['x'][:,0],0).astype(int)
    y_0 = np.round(pos['y'][:,0],0).astype(int)
    x_45 = np.round(pos['x'][:,45],0).astype(int)
    y_45 = np.round(pos['y'][:,45],0).astype(int)
    
    x_0[mask[x_45, y_45] == 0] = -1 
    y_0[mask[x_45, y_45] == 0] = -1 
    x_45[mask[x_45, y_45] == 0] = -1 
    y_45[mask[x_45, y_45] == 0] = -1 
     
    x_0 = x_0[x_0 > -1] 
    y_0 = y_0[y_0 > -1] 
    x_45 = x_45[x_45 > -1] 
    y_45 = y_45[y_45 > -1] 
    
    #Calculate the displacement based on original and contracted images
    analytical_disp = {}
    analytical_disp['x'] = analytical_x[x_45, y_45]
    analytical_disp['y'] = analytical_y[x_45, y_45]
    
    elastix_disp = {}
    elastix_disp['x'] = elastix_x[x_45, y_45]
    elastix_disp['y'] = elastix_y[x_45, y_45]
    
    #linear regression
    lres_x = linregress(analytical_disp['x'],  elastix_disp['x'])
    slope_x = lres_x.slope
    intercept_x = lres_x.intercept
    
    lres_y = linregress(analytical_disp['y'],  elastix_disp['y'])
    slope_y = lres_y.slope
    intercept_y = lres_y.intercept
    
    #r^2 values
    r2_x = linregress(analytical_disp['x'],  elastix_disp['x']).rvalue**2
    r2_y = linregress(analytical_disp['y'],  elastix_disp['y']).rvalue**2
    
    #error
    disp_analytical = np.vstack([analytical_disp['x'], analytical_disp['y']]).T
    disp_elastix = np.vstack([elastix_disp['x'], elastix_disp['y']]).T
    error_elastix = np.linalg.norm(disp_analytical-disp_elastix)
    
    #Display images
    fig, axs = plt.subplots(1, 2)
    axs[0].scatter(analytical_disp['x'], elastix_disp['x'], marker = '.', color='gold')
    x = [-7,10]
    axs[0].plot(x, x, color = 'blue')
    axs[0].set_title('x direction')
    axs[0].annotate('R2=%2.3f' % r2_x, (0.05,0.9), xycoords='axes fraction')
    axs[0].annotate('y=%2.3f' % slope_x + 'x + %2.3f' % intercept_x, (0.05,0.85), xycoords='axes fraction')
    axs[0].set_xlabel("Analytical")
    axs[0].set_ylabel("Elastix")
    
    axs[1].scatter(analytical_disp['y'], elastix_disp['y'], marker = '.', color='gold')
    x2 = [-6,2]
    axs[1].plot(x2, x2, color = 'blue')
    axs[1].set_title('y direction')
    axs[1].annotate('R2=%2.3f' % r2_y, (0.05,0.9), xycoords='axes fraction')
    axs[1].annotate('y=%2.3f' % slope_y + 'x + %2.3f' % intercept_y, (0.05,0.85), xycoords='axes fraction')
    axs[1].set_xlabel("Analytical")
    
    fig.suptitle('Insilico Elastix v Analytical Flow Comparison Plots', fontsize=12)
    
    if save:
        plt.savefig('comparisionPlots/Analytical_v_Elastics_plot.png', dpi = 200)









