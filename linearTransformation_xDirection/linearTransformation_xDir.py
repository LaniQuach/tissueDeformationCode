# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:14:28 2023

@author: laniq
"""

from ITK_TransformImage import displacement_field_elastix, display_save_displacement
from image_to_array import imageToArray, display_bothImages
from transformTestImage import transformImage
from compareImages import compareResults
from strainOperations import secondDeriv, calculate_defGradient, calculate_strain, display_strain, plot_deriv
import os

fldr = 'data/'
out_fldr = fldr + 'out/'
tracking_fldr = fldr + 'contract_analysis/results/'
if not os.path.exists(out_fldr): os.mkdir(out_fldr)

parameterFileName = 'data/parameters_BSpline.txt'
targetImageName = 'data/targetImage.png'
stretchedImageName = 'data/stretchedImage.png'

transformImage(targetImageName, stretchedImageName, fldr)
imageArrays = imageToArray(targetImageName, stretchedImageName)

display_bothImages(imageArrays[0], imageArrays[1], out_fldr + 'originalImages', True)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, True, out_fldr)
display_save_displacement(displacementField, out_fldr + 'dispField', True)

# compare results
compareResults(targetImageName, stretchedImageName, out_fldr + 'elastixResults.npy', out_fldr)

#STRAINS
derivatives = secondDeriv(displacementField[:,:,0], displacementField[:,:,1]*-1)
plot_deriv(derivatives, out_fldr, True)
defGradient = calculate_defGradient(derivatives)
strains = calculate_strain(defGradient)
display_strain(strains, out_fldr, True)




