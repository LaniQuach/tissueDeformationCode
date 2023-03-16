# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:14:28 2023

@author: laniq
"""

from ITK_TransformImage import *
from image_to_array import *
from comparsionPlot_elastix import *
from extractImages import *


mask1 = 'transformingTestImage/Mask_0.png'
mask2 = 'transformingTestImage/Mask_41.png'

analytical_x = np.load('comparisionPlots/analytical_dispx_41.npy')
analytical_y = np.load('comparisionPlots/analytical_dispy_41.npy')

analytical_disp = {}
analytical_disp['x'] = analytical_x[:,:]
analytical_disp['y'] = analytical_y[:,:]

analytical_disp['y'][mask2==0] = np.nan
analytical_disp['x'][mask2==0] = np.nan

parameterFileName = 'data/parameters_BSpline.txt'

# imageArrays = imageToArray(imageFileName1, imageFileName2, mask1, mask2)
imageArrays = [np.load('transformingTestImage/originalArray.npy'),  np.load('transformingTestImage/transformedArray.npy')]
# imageArrays[0] = 
# imageArrays[1] = np.load('transformingTestImage/transformedArray.npy')

display_bothImages(imageArrays[0], imageArrays[1], True)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, mask2)
display_save_displacement(displacementField, True)

#Comparisions
displayComparisionPlots(displacementField[:,:,0], displacementField[:,:,1], analytical_disp['x'], analytical_disp['y'], mask2, True)