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

#convert original image to black and white
image = Image.open('newTestImage.png').convert('L')

plt.figure()
plt.imshow(image)
plt.axis('off')
resizeValue = 100

#add 50 pixels to each side to "stretch" image
width, height = image.size
resized = image.resize((width+resizeValue, height))

#crop transformed image
left = resizeValue/2
right = left+width
top = 0
bottom = height
resizedCrop = resized.crop((left, top, right, bottom))

#turn images into arrays
arrayImg = asarray(resizedCrop)
ogArrayImg = asarray(image)
print(ogArrayImg)

#save arrays to file
np.save('originalArray.npy', ogArrayImg)
np.save('transformedArray.npy', arrayImg)



plt.figure()
plt.imshow(arrayImg, vmin=0, vmax=255)
plt.axis('off')
plt.figure()
plt.imshow(ogArrayImg, vmin=0, vmax=255)
plt.axis('off')


#save the resized image
# resizedCrop.save('newTransformedImage.png')