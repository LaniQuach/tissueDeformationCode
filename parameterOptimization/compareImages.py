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

ogArray = np.load('transformingTestImage/originalArray.npy')
moving_array = np.load('transformingTestImage/transformedArray.npy')
transformedArray = np.load('output/elastixResultImage_Array.npy')

print(np.max(ogArray))
print(np.max(transformedArray))

#crop the two arrays to just the M logo
centerpointX = transformedArray.shape[1]//2
centerpointY = transformedArray.shape[0]//2

# newArray_orig = ogArray[centerpointY-90:centerpointY+75,centerpointX-70:centerpointX+80]
# newArray_transform = transformedArray[centerpointY-90:centerpointY+75,centerpointX-70:centerpointX+80]
newArray_orig = ogArray
newArray_transform = transformedArray

#calculate error
rel_error = np.zeros(newArray_transform.shape)
rel_error[newArray_orig != 0] = (newArray_transform[newArray_orig != 0]
                                 - newArray_orig[newArray_orig != 0] )

outArray = np.subtract(newArray_transform, newArray_orig)

plt.figure()
result = plt.imshow(outArray)
plt.axis('off')
plt.colorbar(result)
plt.savefig('output/errorDifference.png')

avg = np.mean(np.abs(rel_error))
maxVal = np.max(np.abs(outArray))
minVal = np.min(np.abs(outArray))
std = np.std(rel_error)

print("avg error: ", avg)
print("standard deviation", std)
print("max: ", maxVal)
print("min: ", minVal)

fig = plt.figure(figsize=(8, 5))

fig.add_subplot(2,2,1)
fixed = plt.imshow(ogArray, vmin = 0, vmax = 10280)
plt.axis('off')
plt.title("Relaxed", fontsize = 20)

fig.add_subplot(2,2,2)
moving = plt.imshow(moving_array, vmin = 0, vmax = 10280)
plt.axis('off')
plt.title("Contracted", fontsize = 20)

fig.add_subplot(2,2,3)
result = plt.imshow(transformedArray, vmin = 0, vmax = 10280)
plt.axis('off')
plt.title("Elastix Image", fontsize = 20)

fig.add_subplot(2,2,4)
im = plt.imshow(outArray)
plt.colorbar(im)
plt.axis('off')
plt.title("Difference", fontsize = 20)
plt.savefig('output/difference2.png')

plt.axis('off')





