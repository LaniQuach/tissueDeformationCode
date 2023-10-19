# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:14:28 2023

@author: laniq
"""

from ITK_TransformImage import displacement_field_elastix, display_save_displacement
from createImage import createImage
from image_to_array import imageToArray, display_bothImages
from compareImages import compareResults
from strainOperations import secondDeriv, calculate_defGradient, calculate_strain, display_strain, plot_deriv
import os
import numpy as np

fldr = 'data/'
out_fldr = fldr + 'out/'
if not os.path.exists(out_fldr): os.mkdir(out_fldr)

parameterFileName = 'data/parameters_BSpline.txt'
targetImageName = 'data/targetImage.png'
stretchedImageName = 'data/stretchedImage.png'

analytical_gradient = np.asarray(createImage(targetImageName, stretchedImageName)) #F
imageArrays = imageToArray(targetImageName, stretchedImageName)

# display_bothImages(imageArrays[0], imageArrays[1], out_fldr + 'originalImages', True)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, True, out_fldr) #u
display_save_displacement(displacementField, out_fldr + 'dispField', True)

# compare results

#strains using elastix
derivatives_elastix = secondDeriv(displacementField[:,:,0], displacementField[:,:,1]*-1)
plot_deriv(derivatives_elastix, out_fldr, False)
defGradient_elastix = calculate_defGradient(derivatives_elastix) #gradient u
strains_elastix = calculate_strain(defGradient_elastix)
display_strain(strains_elastix, out_fldr + 'elastix_', True)

#strains using analytical
analytical_def = {'x_dx': analytical_gradient[:,:,0,0]-1, 'x_dy': analytical_gradient[:,:,0,1], 'y_dx': analytical_gradient[:,:,1,0], 'y_dy': analytical_gradient[:,:,1,1]-1}
plot_deriv(analytical_def, out_fldr, False)
strains_analyt = calculate_strain(analytical_gradient)
display_strain(strains_analyt, out_fldr + 'analytical_', True)

compareResults(targetImageName, out_fldr + 'analytical_' + 'strain_image.png' , out_fldr + 'elasitx_' + 'strain_image.png', out_fldr)



