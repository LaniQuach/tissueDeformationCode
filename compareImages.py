# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:09:07 2022

@author: Lani
"""

from PIL import Image
import numpy as np
from numpy import asarray
from matplotlib import pyplot as plt
from sklearn import preprocessing

transformedImage = Image.open('output/finalMichLogoUsingBSpline.png')
ogImage = Image.open('output/newTransformedImage.png')

#crop images
# width, height = ogImage.size
# left = (width/6)+10
# right = width - (width/6)
# top = 0
# bottom = height
# ogImage1 = ogImage.crop((left, top, right, bottom))
# ogImage1 = ogImage1.convert('L')

# transformedImage1 = transformedImage.crop((left, top, right, bottom))
# transformedImage1 = transformedImage1.convert('L')


ogArray = asarray(ogImage)
transformedArray = asarray(transformedImage)

# ogArray = ogArray / np.linalg.norm(ogArray)
# transformedArray = transformedArray / np.linalg.norm(transformedArray)
ogArray = preprocessing.normalize(ogArray[:,:,0]) + preprocessing.normalize(ogArray[:,:,1]) + preprocessing.normalize(ogArray[:,:,2])
print(ogArray)
# print(ogArray)
# print(transformedArray)
print(np.max(ogArray))
# outArray = num.subtract(transformedArray, ogArray)
# print(outArray)

# file = open('output/differences.txt', 'w')
# file.write(" ".join(str(x) for x in outArray))
# file.close()



