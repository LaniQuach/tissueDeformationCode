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
import copy

imageFileName1 = 'transformingTestImage/warped_0.png'
imageFileName2 = 'transformingTestImage/warped_41.png'
mask1 = 'transformingTestImage/Mask_0_withPosts.png'
mask2 = 'transformingTestImage/Mask_45_withPosts.png'

imageArrays = img.imageToArray(imageFileName1, imageFileName2, mask1, mask2)

fixed_array = imageArrays[0].astype(np.float32)
moving_array = imageArrays[1].astype(np.float32)

fixed_image = itk.GetImageFromArray(fixed_array)
moving_image = itk.GetImageFromArray(moving_array)

#Initialize values from default parameter list
parameters = {}
parameters['numHistogramBins'] = 32;
parameters['finalGridSpacings'] = 16;
parameters['numResolutions'] = 4;
parameters['numSpatialSamples'] = 2048;
parameters['bSplineInterpOrder'] = 1;
parameters['finalBSplineInterpOrder'] = 3;
parameters['maxNumIterations'] = 500;

parameterFileName = 'data/parameters_BSpline.txt'


resultDispFields = []
for histBins in allParameters.numHistogramBins:
    for gridSpacing in allParameters.finalGridSpacings:
        for numRes in allParameters.numResolutions:
            for numSpatial in allParameters.numSpatialSamples:
                for interpOrder in allParameters.bSplineInterpOrder:
                    for finInterpOrder in allParameters.finalBSplineInterpOrder:
                        try:
                            parameter_object = itk.ParameterObject.New()
            
                            # for testing the correct output
                            parameter_object.AddParameterFile('data/parameters_BSpline.txt')
                            # displacementField = itkTI.displacement_field_elastix(imageArrays[0], imageArrays[1], parameter_object, mask2, True)
                            ###########
            
                            parameter_object.SetParameter("NumberOfHistogramBins", '%i' %histBins)
                            parameters['numHistogramBins'] = histBins
                            
                            parameter_object.SetParameter("FinalGridSpacingInPhysicalUnits", '%i' %gridSpacing)
                            parameters['finalGridSpacings'] = gridSpacing
                            
                            parameter_object.SetParameter("NumberOfResolutions", '%i' %numRes)
                            parameters['numResolutions'] = numRes
                            
                            parameter_object.SetParameter("NumberOfSpatialSamples", '%i' %numSpatial)
                            parameters['numSpatialSamples'] = numSpatial
                            
                            parameter_object.SetParameter("BSplineInterpolationOrder", '%i' %interpOrder)
                            parameters['bSplineInterpOrder'] = interpOrder
                            
                            parameter_object.SetParameter("FinalBSplineInterpolationOrder", '%i' %finInterpOrder)
                            parameters['finalBSplineInterpOrder'] = finInterpOrder
                            
                            displacementField = itkTI.displacement_field_elastix(imageArrays[0], imageArrays[1], parameter_object, mask2, False)
                            # itkTI.display_save_displacement(displacementField, 'outputDispField_parameter', False)
                            
                            data = copy.deepcopy([parameters, displacementField])
                            resultDispFields.append(data)
                            print("works")
                        except:
                            print("parameters of ", histBins, gridSpacing, numRes, numSpatial, interpOrder, finInterpOrder, " don't work")
                            continue
        
np.save('output/allParameter_Results.npy', asarray(resultDispFields, dtype=object))

results = np.load('output/allParameter_Results.npy', allow_pickle=True)
if results.size > 0:
    print(results[:][0])
else:
    print('no values work')



    
    