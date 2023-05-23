#!/usr/bin/env python
# coding: utf-8

# ## 4. Image Registration with initial transform and/or multiple threads

# In this notebook 2 other options of the elastix algorithm are shown: initial transformation and multithreading.
# They're shown together just to reduce the number of example notebooks and 
# thus can be used independently as well as in combination with whichever other functionality
# of the elastix algorithm. 
# 
# Initial transforms are transformations that are done on the moving image before the registration is started.
# 
# Multithreading spreaks for itself and can be used in similar fashion in the transformix algorithm.
# 
# 

# ### Registration

import itk
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from skimage.io import imread

from numpy import asarray

from mpl_toolkits.axes_grid1 import make_axes_locatable

def elastix_transformation(originalArray, movingArray, parameterFileName, saveImageArray, outfldr):
    """
    elastix_transformation utalizes the elastix functionality to attempt to transform 
    moving images back to it's original image

    :param originalArray: original array of the image before contraction/movement
    :param movingArray: array of image during contraction/movement
    :param parameterFileName: name of the file with the ideal parameters for elastix transformation
    :return: the transform parameters of the image that returns the image to it's original state
    """ 
    fixed_array = originalArray.astype(np.float32)
    moving_array = movingArray.astype(np.float32)

    fixed_image = itk.GetImageFromArray(fixed_array)
    moving_image = itk.GetImageFromArray(moving_array)
 
    # Import Default Parameter Map
    parameter_object = itk.ParameterObject.New()
    
    # Load custom parameter maps from .txt file
    parameter_object.AddParameterFile(parameterFileName)
    
    # Load Elastix Image Filter Object
    elastix_object = itk.ElastixRegistrationMethod.New(fixed_image, moving_image)
    elastix_object.SetParameterObject(parameter_object)
    
    # Set additional options
    elastix_object.SetLogToConsole(False)
    
    # Update filter object (required)
    elastix_object.UpdateLargestPossibleRegion()
    
    # Results of Registration
    result_image = elastix_object.GetOutput()
    result_array = itk.GetArrayFromImage(result_image)
    
    if saveImageArray:
        np.save(outfldr + 'elastixResults.npy', result_array)
        
    return elastix_object.GetTransformParameterObject()

def display_save_Image(image, vmin1, vmax1, save):
    """
    display_save_Image displays whatever image you want

    :param image: the image itself
    :param vmin1, vmax1: min and max values of the image to show
    :param save: true or false to save the image
    """ 
    
    plt.figure()
    plt.imshow(image, vmin = vmin1, vmax = vmax1)
    plt.axis('off')
    plt.title(image, 'results')
    if save:
        plt.savefig('output/', image, '.png')
    
def displacement_field_elastix(originalArray, movingArray, parameterFileName, saveImage, outfldr, maskFileName = None):
    """
    displacement_field_elastix creates a displacement field of the two images brought in

    :param originalArray: original array of the image before contraction/movement
    :param movingArray: array of image during contraction/movement
    :param parameterFileName: name of the file with the ideal parameters for elastix transformation
    :return: the displacement field array
    """ 
    resultParameters = elastix_transformation(originalArray, movingArray, parameterFileName, saveImage, outfldr)
    movingImage = itk.GetImageFromArray(movingArray)
    deformation_field = itk.transformix_deformation_field(movingImage, resultParameters)
    defArray = itk.GetArrayFromImage(deformation_field).astype(float)*0.908
    
    if maskFileName is not None:
        mask_1 = np.asarray(imread(maskFileName, as_gray=True))
        defArray[mask_1==0] = np.nan
    
    return defArray

def display_save_displacement(defArray, name, save):
    #Plot images
    fig, axs = plt.subplots(1, 2, sharey=True, figsize=[30,30])
    im3 = axs[1].imshow(defArray[:,:,0], cmap='BrBG')
    
    divider1 = make_axes_locatable(axs[1])
    cax = divider1.new_vertical(size='5%', pad=0.6, pack_start = True)
    fig.add_axes(cax)
    cbar = fig.colorbar(im3, cax = cax, orientation = 'horizontal')
    cbar.set_label('displacement (pixels)', fontsize = 25)
    cbar.ax.tick_params(labelsize=20)

    
    im2 = axs[0].imshow(defArray[:,:,1]*-1, cmap = 'RdYlBu')

    divider = make_axes_locatable(axs[0])
    cax = divider.new_vertical(size='5%', pad=0.6, pack_start = True)
    fig.add_axes(cax)
    cbar2 = fig.colorbar(im2, cax = cax, orientation = 'horizontal')
    cbar2.set_label('displacement (pixels)', fontsize = 25)
    cbar2.ax.tick_params(labelsize=20)
    
    axs[0].axis('off')
    axs[1].axis('off')
    axs[0].set_title('Displacement Field Y', fontsize=30)
    axs[1].set_title('Displacement Field X', fontsize=30)
    
    
    if save: 
        plt.savefig(name + '.png', dpi = 200)

