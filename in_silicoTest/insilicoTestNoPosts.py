# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:14:28 2023

@author: laniq
"""
import numpy as np
from skimage import io
import os

from ITK_TransformImage import displacement_field_elastix, display_save_displacement
from image_to_array import imageToArray, display_bothImages
from comparsionPlot_elastix import displayComparisionPlots


fldr = 'warped_image_v2/'
out_fldr = fldr + 'out/'
tracking_fldr = fldr
if not os.path.exists(out_fldr): os.mkdir(out_fldr)

ts = 45

mask1 = fldr + 'Mask_41.png'
mask2 = fldr + 'Mask_41.png'

parameterFileName = 'data/parameters_BSpline.txt'

image = io.imread(fldr + 'warped.tif')
io.imsave(fldr + 'warped_0.png', image[0])
io.imsave(fldr + 'warped_%i' % ts + '.png', image[ts])

imageArrays = imageToArray(fldr + 'warped_0.png', fldr + 'warped_%i' % ts + '.png')

display_bothImages(imageArrays[0], imageArrays[1], fldr + 'noPosts_frame%i' %ts , True)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, mask2)
display_save_displacement(displacementField, out_fldr + 'dispField_frame%i' % ts + '_noPost', True)
np.save(out_fldr + 'displacement_%i' %ts + '.npy', displacementField)

#Comparisions
displayComparisionPlots(out_fldr, tracking_fldr, fldr, out_fldr, mask2, ts, True)