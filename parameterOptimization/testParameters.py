# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:15:56 2023

@author: laniq
"""

import allParameters
import image_to_array as img
import ITK_TransformImage as itkTI
import itk
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

imageArrays = img.imageToArray(imageFileName1, imageFileName2, mask1, mask2)

fixed_array = imageArrays[0].astype(np.float32)
moving_array = imageArrays[1].astype(np.float32)

fixed_image = itk.GetImageFromArray(fixed_array)
moving_image = itk.GetImageFromArray(moving_array)

#Initialize values from default parameter list
parameters = {}
parameters['finalGridSpacings'] = 16;
parameters['numHistogramBins'] = 32;
parameters['numResolutions'] = 4;
parameters['numSpatialSamples'] = 2048;
parameters['bSplineInterpOrder'] = 1;
parameters['finalBSplineInterpOrder'] = 3;

resultDispFields = []

for histBins in allParameters.numHistogramBins:
    parameter_object = itk.ParameterObject.New()
    parameter_map_bspline = parameter_object.GetDefaultParameterMap("bspline")
    
    parameter_map_bspline['NumberOfHistogramBins'] = str(histBins)
    parameters['numHistogramBins'] = histBins
    
    displacementField = itkTI.displacement_field_elastix(imageArrays[0], imageArrays[1], parameter_map_bspline, mask2)
    itkTI.display_save_displacement(displacementField, 'outputDispField_parameter', False)
    
    data = [parameters, displacementField]
    resultDispFields.append(data)
    
    # plt.figure()
    # plt.imshow(result_image, vmin = 0, vmax = 255)
    # plt.axis('off')
    # display_save_Image(result_image, 0, 255, False)
    
np.save('output/allParameter_Results.npy', asarray(resultDispFields))

    
    