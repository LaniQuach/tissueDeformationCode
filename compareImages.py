# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:09:07 2022

@author: Lani
"""

from PIL import Image
import numpy as num
from numpy import asarray
from matplotlib import pyplot as plt


transformedImage = Image.open('output/finalMichLogoUsingBSpline.png')
ogImage = Image.open('output/newTransformedImage.png')

#crop images
width, height = ogImage.size
left = (width/6)+10
right = width - (width/6)
top = 0
bottom = height
ogImage1 = ogImage.crop((left, top, right, bottom))
ogImage1 = ogImage1.convert('L')

transformedImage1 = transformedImage.crop((left, top, right, bottom))
transformedImage1 = transformedImage1.convert('L')


ogArray = asarray(ogImage1)
transformedArray = asarray(transformedImage1)

outArray = num.subtract(transformedArray, ogArray)

print(outArray)

# file = open('output/differences.txt', 'w')
# file.write(" ".join(str(x) for x in outArray))
# file.close()



