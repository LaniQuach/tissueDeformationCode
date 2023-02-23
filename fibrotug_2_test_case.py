# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 15:20:47 2023

@author: laniq
"""

####TEST CASE####
from ITK_TransformImage import *
from image_to_array import *
from comparsionPlot import *

imageFileName1 = 'fibrotug_2/transformingTestImage/A_0.41_01-14-2_rc_ds_0.png'
imageFileName2 = 'fibrotug_2/transformingTestImage/A_0.41_01-14-2_rc_ds_25.png'
mask1 = 'fibrotug_2/transformingTestImage/Mask_0.png'
mask2 = 'fibrotug_2/transformingTestImage/Mask_25.png'

parameterFileName = 'fibrotug_2/data/parameters_BSpline.txt'

imageArrays = imageToArray(imageFileName1, imageFileName2, mask1, mask2)
display_bothImages(imageArrays[0], imageArrays[1], True)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, mask2)
display_save_displacement(displacementField)

#Comparisions
displayComparisionPlots(displacementField[:,:,0], displacementField[:,:,1], mask2, True)
