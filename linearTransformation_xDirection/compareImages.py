# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:09:07 2022

@author: Lani
"""

from PIL import Image
import numpy as np
from numpy import asarray
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
import itk

#read in all images
fixed_image = itk.imread('output/originalImage.png', itk.F)
moving_image = itk.imread('output/stretchedImage.png', itk.F)
result_image = itk.imread('output/elastixResults.png', itk.F)

#transform to arrays
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
newArray_orig = ogArray[40:250,100:345]
newArray_transform = transformedArray[40:250,100:345]

#calculate error
rel_error = np.zeros(newArray_transform.shape)
rel_error[newArray_orig != 0] = (newArray_transform[newArray_orig != 0]
                                 - newArray_orig[newArray_orig != 0] )

outArray = np.subtract(newArray_transform, newArray_orig)

plt.figure()
result = plt.imshow(outArray)
plt.colorbar(result)
plt.axis('off')
plt.savefig('output/errorDifference.png')

avg = np.mean(rel_error)
maxVal = np.max(np.abs(outArray))
minVal = np.min(np.abs(outArray))
std = np.std(rel_error)

print("avg error: ", avg)
print("standard deviation", std)
print("max: ", maxVal)
print("min: ", minVal)

# fig = plt.figure(figsize=(9, 7))

# fig.add_subplot(2,2,1)
# fixed = plt.imshow(fixed_image)
# plt.colorbar(fixed)
# plt.axis('off')
# plt.title("Original Image")

# fig.add_subplot(2,2,2)
# moving = plt.imshow(moving_image)
# plt.colorbar(moving)
# plt.axis('off')
# plt.title("Stretched Image")

# fig.add_subplot(2,2,3)
# result = plt.imshow(transformedArray)
# plt.colorbar(result)
# plt.axis('off')
# plt.title("Elastix Image")

# fig.add_subplot(2,2,4)
# im = plt.imshow(outArray)
# plt.colorbar(im)
# plt.axis('off')
# plt.title("Difference")
# plt.savefig('output/difference2.png')

# plt.axis('off')





