#!/usr/bin/env python
# coding: utf-8

# ### Registration
import itk
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from skimage.io import imread
from numpy import asarray
from mpl_toolkits.axes_grid1 import make_axes_locatable

#load both arrays
def elastix_transformation(originalArray, movingArray, parameterMap):
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
    parameter_object.AddParameterMap(parameterMap)
    
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
    plt.title(image, ' results')
    if save:
        plt.savefig('fibro_tug/output/', image, '.png')
    
def displacement_field_elastix(originalArray, movingArray, parameterMap, maskFileName):
    """
    displacement_field_elastix creates a displacement field of the two images brought in

    :param originalArray: original array of the image before contraction/movement
    :param movingArray: array of image during contraction/movement
    :param parameterFileName: name of the file with the ideal parameters for elastix transformation
    :return: the displacement field array
    """ 
    resultParameters = elastix_transformation(originalArray, movingArray, parameterMap)
    movingImage = itk.GetImageFromArray(movingArray)
    deformation_field = itk.transformix_deformation_field(movingImage, resultParameters)
    defArray = itk.GetArrayFromImage(deformation_field).astype(float)*0.908
    
    mask_1 = np.asarray(imread(maskFileName, as_gray=True))
    defArray[mask_1==0] = np.nan
    
    return defArray

def displacement_field_elastix_withoutmask(originalArray, movingArray, parameterFileName):
    """
    displacement_field_elastix creates a displacement field of the two images brought in

    :param originalArray: original array of the image before contraction/movement
    :param movingArray: array of image during contraction/movement
    :param parameterFileName: name of the file with the ideal parameters for elastix transformation
    :return: the displacement field array
    """ 
    resultParameters = elastix_transformation(originalArray, movingArray, parameterFileName)
    movingImage = itk.GetImageFromArray(movingArray)
    deformation_field = itk.transformix_deformation_field(movingImage, resultParameters)
    defArray = itk.GetArrayFromImage(deformation_field).astype(float)*0.908
    
    return defArray

def display_save_displacement(defArray, name, save):
    #Plot images
    fig, axs = plt.subplots(1, 2, sharey=True, figsize=[30,30])
    im3 = axs[1].imshow(defArray[:,:,0], vmin = -10, vmax = 10)
    #, vmin = -15, vmax = 15
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im3, cax=cax, orientation='vertical');
    cbar.set_label('displacement (pixels)', fontsize = 25)
    cbar.ax.tick_params(labelsize=30)
    
    
    im2 = axs[0].imshow(defArray[:,:,1]*-1, vmin = -7, vmax = 7)
    #, vmin = -5, vmax = 5
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar2 = fig.colorbar(im2, cax=cax, orientation='vertical');
    cbar2.ax.tick_params(labelsize=30)
    cbar2.set_label('displacement (pixels)', fontsize = 25)
    axs[0].axis('off')
    axs[1].axis('off')
    axs[0].set_title('Displacement Field Y', fontsize=30)
    axs[1].set_title('Displacement Field X', fontsize=30)
    
    if save: 
        plt.savefig('output/' + name + '.png', dpi = 200)
    



