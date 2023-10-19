# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:14:28 2023

@author: laniq
"""

from ITK_TransformImage import displacement_field_elastix, display_save_displacement
from createImage import createImage
from image_to_array import imageToArray, display_bothImages
from compareImages import compareStrain
from strainOperations import secondDeriv, calculate_defGradient, calculate_strain, display_strain, plot_deriv
import os
import numpy as np
from matplotlib import pyplot as plt

fldr = 'data/'
out_fldr = fldr + 'out/'
if not os.path.exists(out_fldr): os.mkdir(out_fldr)

parameterFileName = 'data/parameters_BSpline.txt'
targetImageName = 'data/targetImage.png'
stretchedImageName = 'data/stretchedImage.png'

vmin = [-0.12, -0.3, -0.3, -0.05]
vmax = [0.12, 0.3, 0.3, 0.15]

analytical_gradient, I_disp, J_disp, E = createImage(targetImageName, stretchedImageName) #F
analytical_gradient = np.asarray(analytical_gradient)
imageArrays = imageToArray(targetImageName, stretchedImageName)

# # display_bothImages(imageArrays[0], imageArrays[1], out_fldr + 'originalImages', True)
display_strain(E, out_fldr, vmin, vmax, False)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, True, out_fldr) #u
display_save_displacement(displacementField, out_fldr + 'dispField', True)

#strains using elastix
derivatives_elastix = secondDeriv(displacementField[:,:,1], displacementField[:,:,0])
# plot_deriv(derivatives_elastix, out_fldr, False)
defGradient_elastix = calculate_defGradient(derivatives_elastix) #gradient u
defGradient_elastix = np.round(defGradient_elastix, decimals = 2)

#elastix strain
strains_elastix = calculate_strain(defGradient_elastix)
display_strain(strains_elastix, out_fldr + 'elastix_', vmin, vmax, True)

#strains using analytical finite difference
analytical_finiteDiff = secondDeriv(I_disp, J_disp)
# plot_deriv(analytical_finiteDiff, out_fldr + 'analytical_', False)
defGradient_analytical = calculate_defGradient(analytical_finiteDiff) #gradient u
defGradient_analytical = np.round(defGradient_analytical, decimals = 2)
strains_analytical_finiteDiff = calculate_strain(defGradient_analytical)
display_strain(strains_analytical_finiteDiff, out_fldr + 'elastix_', vmin, vmax, True)

#strains using analytical exact
analytical_def = {'x_dx': analytical_gradient[:,:,0,0]-1, 'x_dy': analytical_gradient[:,:,0,1], 'y_dx': analytical_gradient[:,:,1,0], 'y_dy': analytical_gradient[:,:,1,1]-1}
# plot_deriv(analytical_def, out_fldr, False)
strains_analyt = calculate_strain(analytical_gradient)
display_strain(strains_analyt, out_fldr + 'analytical_', vmin, vmax, True)

compareStrain(strains_elastix, strains_analyt, vmin, vmax, out_fldr)



