# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:09:07 2022

@author: Lani
"""

from PIL import Image
import numpy as np
from numpy import asarray
from matplotlib import pyplot as plt
import itk

fixed_image = itk.imread('transformingTestImage/newTestImage.png', itk.F)
result_image = itk.imread('output/result_image_11_7.nii', itk.F)
# transformedImage = Image.open('output/finalMichLogoUsingBSpline.png')
# ogImage = Image.open('output/newTransformedImage.png')

# #crop images
# width, height = ogImage.size
# left = (width/9)
# right = width - (width/9)
# top = 0
# bottom = height

# ogImage1 = ogImage.crop((left, top, right, bottom))
# transformedImage1 = transformedImage.crop((left, top, right, bottom))

# ogArray = np.array(ogImage1)
# transformedArray = np.array(transformedImage1)

ogArray = itk.GetArrayFromImage(fixed_image)
transformedArray = itk.GetArrayFromImage(result_image)

#Normalize values to 0-1
maxValOG = np.max(ogArray)
minValOG = np.min(ogArray)
ogArray = (ogArray-minValOG)/(maxValOG-minValOG)  

maxVal_Transformed = np.max(transformedArray)
minVal_Transformed = np.min(transformedArray)
transformedArray = (transformedArray-minVal_Transformed)/(maxVal_Transformed-minVal_Transformed)  

newArray_orig = ogArray[:,40:290]
newArray_transform = transformedArray[:,40:290]

outArray = np.subtract(newArray_transform, newArray_orig)
print(outArray)
im = plt.imshow(outArray)
plt.colorbar(im)

file = open('output/test.txt', 'w')
file.write(" ".join(str(x) for x in outArray))
file.close()



