# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:14:28 2023

@author: laniq
"""

import ITK_TransformImage as itk
import image_to_array as img
import comparsionPlot_elastix as cmp
import extractImages as exti

tissueName = 'transformingTestImage/warped'

contractedFrame = 45
exti.getIndividualImages(tissueName, 0)
exti.getIndividualImages(tissueName, contractedFrame)

imageFileName1 = tissueName + '_%i' % 0 + '.png'
imageFileName2 = tissueName + '_%i' % contractedFrame + '.png'

mask1 = 'transformingTestImage/Mask_0_withPosts.png'
mask2 = 'transformingTestImage/Mask_45_withPosts.png'

analytical_x = cmp.np.load('comparisionPlots' + '/analytical_dispx_%i' % contractedFrame + '_withPosts.npy')
analytical_y = cmp.np.load('comparisionPlots' + '/analytical_dispy_%i' % contractedFrame + '_withPosts.npy')

analytical_disp = {}
analytical_disp['x'] = analytical_x[:,:]
analytical_disp['y'] = analytical_y[:,:]

analytical_disp['y'][mask2==0] = cmp.np.nan
analytical_disp['x'][mask2==0] = cmp.np.nan

parameterFileName = 'data/parameters_BSpline.txt'

imageArrays = img.imageToArray(imageFileName1, imageFileName2, mask1, mask2)
img.display_bothImages(imageArrays[0], imageArrays[1], 'originalContractedWithPosts', True)

#Elastix functionality
displacementField = itk.displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, mask2)
itk.display_save_displacement(displacementField, 'outputDispField_withPosts', True)

#Comparisions
cmp.displayComparisionPlots(displacementField[:,:,0], displacementField[:,:,1], analytical_disp['x'], analytical_disp['y'], mask2, contractedFrame, True)