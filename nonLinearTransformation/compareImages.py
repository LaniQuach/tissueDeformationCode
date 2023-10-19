# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:09:07 2022

@author: Lani
"""

from PIL import Image
import numpy as np
from numpy import asarray
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
import itk

def compareResults(targetImageName, stretchedImageName, elastixArrayName, outfldr):
    #read in all images
    fixed_image = itk.imread(targetImageName, itk.F)
    moving_image = itk.imread(stretchedImageName, itk.F)
    result_image = itk.imread(elastixArrayName, itk.F)
    
    #transform to arrays
    ogArray = itk.GetArrayFromImage(fixed_image)
    transformedArray = itk.GetArrayFromImage(result_image)
    movingArray = itk.GetArrayFromImage(moving_image)
    
    #Normalize values to 0-1
    maxValOG = np.max(ogArray)
    minValOG = np.min(ogArray)
    ogArray = (ogArray-minValOG)/(maxValOG-minValOG)  
    
    maxVal_Transformed = np.max(transformedArray)
    minVal_Transformed = np.min(transformedArray)
    transformedArray = (transformedArray-minVal_Transformed)/(maxVal_Transformed-minVal_Transformed)  
    #####
    
    #crop the two arrays to just the M logo
    centerpointX = transformedArray.shape[1]//2
    centerpointY = transformedArray.shape[0]//2

    newArray_orig = ogArray[centerpointY-90:centerpointY+75,centerpointX-110:centerpointX+120]
    newArray_transform = transformedArray[centerpointY-90:centerpointY+75,centerpointX-110:centerpointX+120]

    # newArray_orig = ogArray
    # newArray_transform = transformedArray

    #calculate error
    rel_error = np.zeros(newArray_transform.shape)
    rel_error[newArray_orig != 0] = (newArray_transform[newArray_orig != 0]
                                     - newArray_orig[newArray_orig != 0] )

    outArray = np.subtract(newArray_transform, newArray_orig)

    plt.figure()
    result = plt.imshow(outArray, cmap='plasma')
    plt.axis('off')
    cbar = plt.colorbar(result)
    cbar.set_label(label = 'error', size = 15)
    plt.title('Relative Error')
    plt.savefig(outfldr + 'errorDifference.png', dpi = 200)

    avg = np.mean(np.abs(rel_error))
    maxVal = np.max(np.abs(outArray))
    minVal = np.min(np.abs(outArray))
    std = np.std(rel_error)

    print("avg error: ", avg)
    print("standard deviation", std)
    print("max: ", maxVal)
    print("min: ", minVal)

    fig = plt.figure(figsize=(24, 18))

    fig.add_subplot(2,2,1)
    plt.imshow(ogArray, vmin = 0, vmax = 1)
    plt.axis('off')
    plt.title("Original Image")

    fig.add_subplot(2,2,2)
    plt.imshow(movingArray, vmin = 0, vmax = 255)
    plt.axis('off')
    plt.title("Stretched Image")

    fig.add_subplot(2,2,3)
    result = plt.imshow(transformedArray, vmin = 0, vmax = 1)
    plt.axis('off')
    plt.title("Elastix Image")

    fig.add_subplot(2,2,4)
    im = plt.imshow(outArray, cmap='plasma')
    plt.colorbar(im)
    plt.axis('off')
    plt.title("Difference")
    plt.savefig(outfldr + 'allResults.png')

    plt.axis('off')
    
def compareStrain(elastix_strain, analytical_strain, vmin, vmax, outfldr):
    #Normalize values to 0-1
    analytical_strain[:,:,0,0] = (analytical_strain[:,:,0,0]-vmin[0])/(vmax[0]-vmin[0])  
    analytical_strain[:,:,0,1] = (analytical_strain[:,:,0,1]-vmin[1])/(vmax[1]-vmin[1])  
    analytical_strain[:,:,1,0] = (analytical_strain[:,:,1,0]-vmin[2])/(vmax[2]-vmin[2])  
    analytical_strain[:,:,1,1] = (analytical_strain[:,:,1,1]-vmin[3])/(vmax[3]-vmin[3])  

    elastix_strain[:,:,0,0] = (elastix_strain[:,:,0,0]-vmin[0])/(vmax[0]-vmin[0])  
    elastix_strain[:,:,0,1] = (elastix_strain[:,:,0,1]-vmin[1])/(vmax[1]-vmin[1])  
    elastix_strain[:,:,1,0] = (elastix_strain[:,:,1,0]-vmin[2])/(vmax[2]-vmin[2])  
    elastix_strain[:,:,1,1] = (elastix_strain[:,:,1,1]-vmin[3])/(vmax[3]-vmin[3])  
    #####

    #calculate error
    rel_error_xx = np.zeros(elastix_strain[:,:,0,0].shape)
    rel_error_xx = (elastix_strain[:,:,0,0]
                                      - analytical_strain[1:-1,1:-1,0,0])
    
    rel_error_xy = np.zeros(elastix_strain[:,:,0,1].shape)
    rel_error_xy = (elastix_strain[:,:,0,1]
                                      - analytical_strain[1:-1,1:-1,0,1])
    
    rel_error_yx = np.zeros(elastix_strain[:,:,1,0].shape)
    rel_error_yx = (elastix_strain[:,:,1,0]
                                      - analytical_strain[1:-1,1:-1,1,0])
    
    rel_error_yy = np.zeros(elastix_strain[:,:,1,1].shape)
    rel_error_yy = (elastix_strain[:,:,1,1]
                                      - analytical_strain[1:-1,1:-1,1,1])

    avg = [np.mean(np.abs(rel_error_xx)), np.mean(np.abs(rel_error_xy)), np.mean(np.abs(rel_error_yx)), np.mean(np.abs(rel_error_yy))]
    maxVal = [np.max(rel_error_xx), np.max(rel_error_xy), np.max(rel_error_yx), np.max(rel_error_yy)]
    minVal = [np.min(rel_error_xx), np.min(rel_error_xy), np.min(rel_error_yx), np.min(rel_error_yy)]
    std = [np.std(rel_error_xx), np.std(rel_error_xy), np.std(rel_error_yx), np.std(rel_error_yy)]

    print("avg error: ", avg)
    print("standard deviation", std)
    print("max: ", maxVal)
    print("min: ", minVal)

    fig = plt.figure(figsize=(24, 18))

    fig.add_subplot(2,2,1)
    result = plt.imshow(rel_error_xx)
    cbar1 = plt.colorbar(result)
    plt.axis('off')
    plt.title("xx", fontsize=30)

    fig.add_subplot(2,2,2)
    result2 = plt.imshow(rel_error_xy)
    cbar2 = cbar3 = plt.colorbar(result2)
    plt.axis('off')
    plt.title("xy", fontsize=30)

    fig.add_subplot(2,2,3)
    result3 = plt.imshow(rel_error_yx)
    plt.colorbar(result3)
    plt.axis('off')
    plt.title("yx", fontsize=30)

    fig.add_subplot(2,2,4)
    im = plt.imshow(rel_error_yy, cmap='plasma')
    plt.colorbar(im)
    cbar1.ax.tick_params(labelsize=30)

    plt.axis('off')
    plt.title("yy", fontsize=30)
    # plt.savefig(outfldr + 'allResults.png')

    plt.axis('off')






