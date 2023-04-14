# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 10:30:22 2022

@author: Lani
"""

from PIL import Image
from numpy import asarray
import numpy as np
from matplotlib import pyplot as plt

def transformImage(imageName, outputImageName, path):
    #convert original image to black and white
    image = Image.open(imageName).convert('L')
    
    plt.figure()
    # plt.imshow(image)
    # plt.axis('off')
    
    #add 50 pixels to each side to "stretch" image
    resizeValue = 100
    width, height = image.size
    resized = image.resize((width+resizeValue, height))
    
    #crop transformed image
    left = resizeValue/2
    right = left+width
    top = 0
    bottom = height
    resizedCrop = resized.crop((left, top, right, bottom))
    
    resizedCrop.save(outputImageName)

    plt.figure()
    plt.imshow(asarray(image), vmin=0, vmax=255)
    plt.axis('off')
    plt.figure()
    plt.imshow(asarray(resizedCrop), vmin=0, vmax=255)
    plt.axis('off')
    