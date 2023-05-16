# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 14:50:04 2023

@author: laniq
"""
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from scipy.stats import linregress


def displayComparisionPlots(elastix_fldr, tracking_fldr, analytical_fldr, out_fldr, movingMaskFileName, frame, save):
    """
    displayComparisionPlots displays two different sets of displacement data to one another

    :param displacementArray_x, displacementArray_y: first set of displacement data in the x and y fields
    :param movingMaskFileName: file name of the mask to apply to image
    :param save: whether or not to save the resulting displayed image
    """
    beat = 0

    #Emma's data
    pos = {}
    pos['x'] = np.loadtxt(tracking_fldr + 'beat%i' % beat + '_col.txt')
    pos['y'] = np.loadtxt(tracking_fldr + 'beat%i' % beat + '_row.txt')
    tracking_x = pos['x'][:,frame] - pos['x'][:,0]
    tracking_y = pos['y'][:,frame] - pos['y'][:,0]

    #Elastix data
    displacement = np.load(elastix_fldr + 'displacement_%i' %frame + '.npy')
    displacementArray_x = displacement[:,:,0]
    displacementArray_y = displacement[:,:,1]
    elastix_x = displacementArray_x.T
    elastix_y = (displacementArray_y).T

    #Analytical data
    analyticalDisp_x = np.load(analytical_fldr + 'analytical_dispx_%i' %frame + '_withPosts.npy')
    analyticalDisp_y = np.load(analytical_fldr + 'analytical_dispy_%i' %frame + '_withPosts.npy')
    analytical_x = analyticalDisp_x.T
    analytical_y = analyticalDisp_y.T

    mask = np.asarray(imread(movingMaskFileName, as_gray=True)).T

    #use these values to find coordinates to track on elastix's displacement field
    x_0 = np.round(pos['x'][:,0],0).astype(int)
    y_0 = np.round(pos['y'][:,0],0).astype(int)
    x_contracted = np.round(pos['x'][:,frame],0).astype(int)
    y_contracted = np.round(pos['y'][:,frame],0).astype(int)
    mask_track = mask[x_contracted, y_contracted]

    x_0[mask[x_contracted, y_contracted] == 0] = -1
    y_0[mask[x_contracted, y_contracted] == 0] = -1
    x_contracted[mask[x_contracted, y_contracted] == 0] = -1
    y_contracted[mask[x_contracted, y_contracted] == 0] = -1

    x_0 = x_0[x_0 > -1]
    y_0 = y_0[y_0 > -1]
    x_contracted = x_contracted[x_contracted > -1]
    y_contracted = y_contracted[y_contracted > -1]

    #Calculate the displacement based on original and contracted images
    analytical_disp = {}
    analytical_disp['x'] = analytical_x[x_contracted, y_contracted]
    analytical_disp['y'] = analytical_y[x_contracted, y_contracted]

    elastix_disp = {}
    elastix_disp['x'] = elastix_x[x_contracted, y_contracted]
    elastix_disp['y'] = elastix_y[x_contracted, y_contracted]

    tracking_disp = {}
    tracking_disp['x'] = tracking_x[mask_track != 0]
    tracking_disp['y'] = tracking_y[mask_track != 0]

    #linear regression
    lres_x = linregress(analytical_disp['x'],  elastix_disp['x'])
    slope_x = lres_x.slope
    intercept_x = lres_x.intercept

    lres_y = linregress(analytical_disp['y'],  elastix_disp['y'])
    slope_y = lres_y.slope
    intercept_y = lres_y.intercept

    #r^2 values elastix
    r2_x_elastix = linregress(analytical_disp['x'],  elastix_disp['x']).rvalue**2
    r2_y_elastix = linregress(analytical_disp['y'],  elastix_disp['y']).rvalue**2
    
    #r^2 values EMMA
    r2_x_emma = linregress(analytical_disp['x'],  tracking_disp['x']).rvalue**2
    r2_y_emma = linregress(analytical_disp['y'],  tracking_disp['y']).rvalue**2

    #error
    disp_analytical = np.vstack([analytical_disp['x'], analytical_disp['y']]).T
    disp_elastix = np.vstack([elastix_disp['x'], elastix_disp['y']]).T
    error_elastix = np.linalg.norm(disp_analytical-disp_elastix)

    #Display images
    fig, axs = plt.subplots(1, 2)
    axs[0].scatter(analytical_disp['x'], elastix_disp['x'], marker = '.', color='gold')
    axs[0].scatter(analytical_disp['x'], tracking_disp['x'], marker = '.', color='red')
    x = [-7,6]
    axs[0].plot(x, x, color = 'blue')
    axs[0].set_title('x direction')
    axs[0].annotate('R2 (Elastix)=%2.3f' % r2_x_elastix, (0.05,0.9), xycoords='axes fraction')
    axs[0].annotate('R2 (Optical Flow)=%2.3f' % r2_x_emma, (0.05,0.95), xycoords='axes fraction')
    # axs[0].annotate('y=%2.3f' % slope_x + 'x + %2.3f' % intercept_x, (0.05,0.85), xycoords='axes fraction')
    axs[0].legend(["Elastix", "Optical Flow"], loc = 'lower right', frameon=False)

    axs[0].set_xlabel("Analytical")
    axs[0].set_ylabel("Measured")
    
    axs[1].scatter(analytical_disp['y'], elastix_disp['y'], marker = '.', color='gold')
    axs[1].scatter(analytical_disp['y'], tracking_disp['y'], marker = '.', color='red')
    x2 = [-6,2]
    axs[1].plot(x2, x2, color = 'blue')
    axs[1].set_title('y direction')
    axs[1].annotate('R2 (Elastix)=%2.3f' % r2_y_elastix, (0.05,0.9), xycoords='axes fraction')
    axs[1].annotate('R2 (Optical Flow)=%2.3f' % r2_y_emma, (0.05,0.95), xycoords='axes fraction')
    axs[1].legend(["Elastix", "Optical Flow"], loc = 'lower right', frameon=False)
    # axs[1].annotate('y=%2.3f' % slope_y + 'x + %2.3f' % intercept_y, (0.05,0.85), xycoords='axes fraction')
    axs[1].set_xlabel("Analytical")

    # fig.suptitle('Insilico Elastix v Analytical Flow Comparison Plots', fontsize=12)
    
    if save:
        plt.savefig(out_fldr + 'Analytical_v_Elastics_plot_%i' % frame + '.png', dpi = 200)

    
    #error
    # disp_analytical = np.vstack([analytical_disp['x'], analytical_disp['y']]).T
    # disp_elastix = np.vstack([elastix_disp['x'], elastix_disp['y']]).T
    # error_elastix = np.linalg.norm(disp_analytical-disp_elastix)
    






