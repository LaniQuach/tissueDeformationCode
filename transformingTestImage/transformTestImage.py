# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 10:30:22 2022

@author: Lani
"""

from PIL import Image

image = Image.open('newTestImage.png')

resizeValue = 100

width, height = image.size
resized = image.resize((width+resizeValue, height))
image.show()

#crop transformed image
left = resizeValue/2
right = left+width
top = 0
bottom = height
resizedCrop = resized.crop((left, top, right, bottom))
resizedCrop.show()

resizedCrop.save('newTransformedImage.png')