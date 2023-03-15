# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:15:56 2023

@author: laniq
"""

import allParameters
from image_to_array import *
from ITK_TransformImage import display_save_Image, displacement_field_elastix
import itk
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from skimage.io import imread
from numpy import asarray
from mpl_toolkits.axes_grid1 import make_axes_locatable

imageFileName1 = 'transformingTestImage/warped_0.png'
imageFileName2 = 'transformingTestImage/warped_41.png'
mask1 = 'transformingTestImage/Mask_0.png'
mask2 = 'transformingTestImage/Mask_41.png'

imageArrays = imageToArray(imageFileName1, imageFileName2, mask1, mask2)

fixed_array = imageArrays[0].astype(np.float32)
moving_array = imageArrays[1].astype(np.float32)

fixed_image = itk.GetImageFromArray(fixed_array)
moving_image = itk.GetImageFromArray(moving_array)

for metric in allParameters.metrics:
    parameter_object = itk.ParameterObject.New()
    parameter_map_bspline = parameter_object.GetDefaultParameterMap("bspline")
    parameter_map_bspline['Metric'] = [metric]
    parameter_object.AddParameterMap(parameter_map_bspline)
    
    elastix_object = itk.ElastixRegistrationMethod.New(fixed_image, moving_image)
    elastix_object.SetParameterObject(parameter_object)

    elastix_object.SetLogToConsole(False)

    elastix_object.UpdateLargestPossibleRegion()

    # Results of Registration
    result_image = elastix_object.GetOutput()
    result_array = itk.GetArrayFromImage(result_image)
    transformParameter = elastix_object.GetTransformParameterObject()
    
    displacement_field_elastix(movingImage, transformParameter, mask2)
    np.save('output/metrics/metric_' + metric + '.npy', defArray[:,:,0])
    
    # plt.figure()
    # plt.imshow(result_image, vmin = 0, vmax = 255)
    # plt.axis('off')
    # display_save_Image(result_image, 0, 255, False)
    

    
    