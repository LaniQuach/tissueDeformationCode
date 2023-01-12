# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 10:30:22 2022

@author: Lani
"""

from PIL import Image
from numpy import asarray
import numpy as np
from skimage.color import rgb2gray
from matplotlib import pyplot as plt

image = Image.open('newTestImage.png').convert('L')

resizeValue = 100

width, height = image.size
resized = image.resize((width+resizeValue, height))
# image.show()

#crop transformed image
left = resizeValue/2
right = left+width
top = 0
bottom = height
resizedCrop = resized.crop((left, top, right, bottom))
# resizedCrop.show()
arrayImg = asarray(resizedCrop)
ogArrayImg = asarray(image)

np.save('original.npy', ogArrayImg)
np.save('transformed.npy', arrayImg)

plt.figure()
plt.imshow(arrayImg, vmin=0, vmax=255)

plt.figure()
plt.imshow(ogArrayImg, vmin=0, vmax=255)

# resizedCrop.save('newTransformedImage.png')