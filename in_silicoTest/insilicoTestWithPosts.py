# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:14:28 2023

@author: laniq
"""

from ITK_TransformImage import *
from image_to_array import *
from comparsionPlot_elastix import *
from extractImages import *

tissueName = 'transformingTestImage/warped'

contractedFrame = 45
getIndividualImages(tissueName, 0)
getIndividualImages(tissueName, contractedFrame)

imageFileName1 = tissueName + '_%i' % 0 + '.png'
imageFileName2 = tissueName + '_%i' % contractedFrame + '.png'

mask1 = 'transformingTestImage/Mask_0_withPosts.png'
mask2 = 'transformingTestImage/Mask_45_withPosts.png'

analytical_x = np.load('comparisionPlots/analytical_dispx_45_withPosts.npy')
analytical_y = np.load('comparisionPlots/analytical_dispy_45_withPosts.npy')

analytical_disp = {}
analytical_disp['x'] = analytical_x[:,:]
analytical_disp['y'] = analytical_y[:,:]

analytical_disp['y'][mask2==0] = np.nan
analytical_disp['x'][mask2==0] = np.nan

parameterFileName = 'data/parameters_BSpline.txt'

imageArrays = imageToArray(imageFileName1, imageFileName2, mask1, mask2)
display_bothImages(imageArrays[0], imageArrays[1], True)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, mask2)
display_save_displacement(displacementField, True)

#Comparisions
displayComparisionPlots(displacementField[:,:,0], displacementField[:,:,1], analytical_disp['x'], analytical_disp['y'], mask2, True)