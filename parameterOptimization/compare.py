# -*- coding: utf-8 -*-
"""
Created on Tue May 16 09:33:36 2023

@author: laniq
"""
import numpy as np
from skimage import io
import os

from ITK_TransformImage_orig import displacement_field_elastix, display_save_displacement
from image_to_array import imageToArray, display_bothImages
from comparsionPlot import displayComparisionPlots

# fldr = 'warped_image_post_v2/'
# out_fldr = fldr + 'out/'
# tracking_fldr = fldr + 'contract_analysis/results/'
# if not os.path.exists(out_fldr): os.mkdir(out_fldr)

ts = 41

mask1 = 'transformingTestImage/Mask_0_withPosts.png'
mask2 = 'transformingTestImage/Mask_45_withPosts.png'

analytical_x = np.load('comparisionPlots/' + 'analytical_dispx_%i' % ts + '_withPosts.npy').T
analytical_y = np.load('comparisionPlots/' + 'analytical_dispy_%i' % ts + '_withPosts.npy').T

analytical_disp = {}
analytical_disp['x'] = analytical_x[:,:]
analytical_disp['y'] = analytical_y[:,:]

analytical_disp['y'][mask2==0] = np.nan
analytical_disp['x'][mask2==0] = np.nan

parameterFileName = 'data/parameters_BSpline_optimized.txt'

# image = io.imread(fldr + 'warped.tif')
# io.imsave(fldr + 'warped_0_posts.png', image[0])
# io.imsave(fldr + 'warped_%i_posts' % ts + '.png', image[ts])

imageArrays = imageToArray('transformingTestImage/' + 'warped_0.png', 'transformingTestImage/' + 'warped_41.png', mask1, mask2)

display_bothImages(imageArrays[0], imageArrays[1], 'transformingTestImage/' + 'original_images.png' , False)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, mask2)
display_save_displacement(displacementField, 'output/' + 'dispField_frame%i' % ts + '_Post.png', True)
np.save('comparisionPlots/' + 'displacement_%i' %ts + '.npy', displacementField)


#Comparisions
displayComparisionPlots('comparisionPlots/', 'comparisionPlots/', 'comparisionPlots/', 'comparisionPlots/', mask2, ts, True)