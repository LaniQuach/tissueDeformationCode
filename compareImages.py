# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:09:07 2022

@author: Lani
"""

from PIL import Image
import numpy as np
from numpy import asarray
from matplotlib import pyplot as plt
import pandas as pd
import itk

fixed_image = itk.imread('transformingTestImage/newTestImage.png', itk.F)
result_image = itk.imread('output/result_image_11_7.nii', itk.F)

ogArray = itk.GetArrayFromImage(fixed_image)
transformedArray = itk.GetArrayFromImage(result_image)

#Normalize values to 0-1
maxValOG = np.max(ogArray)
minValOG = np.min(ogArray)
ogArray = (ogArray-minValOG)/(maxValOG-minValOG)  

maxVal_Transformed = np.max(transformedArray)
minVal_Transformed = np.min(transformedArray)
transformedArray = (transformedArray-minVal_Transformed)/(maxVal_Transformed-minVal_Transformed)  
#####

#crop the two arrays to just the M logo
newArray_orig = ogArray[10:210,40:290]
newArray_transform = transformedArray[10:210,40:290]

rel_error = np.zeros(newArray_transform.shape)
rel_error[newArray_orig != 0] = (newArray_transform[newArray_orig != 0] - newArray_orig[newArray_orig != 0] )/newArray_orig[newArray_orig != 0]

outArray = np.subtract(newArray_transform, newArray_orig)

avg = np.mean(rel_error)
maxVal = np.max(np.abs(outArray))
minVal = np.min(np.abs(outArray))
std = np.std(rel_error)

print("avg: ", avg)
print("standard deviation", std)
print("max: ", maxVal)
print("min: ", minVal)

im = plt.imshow(outArray)
plt.colorbar(im)
plt.axis('off')





