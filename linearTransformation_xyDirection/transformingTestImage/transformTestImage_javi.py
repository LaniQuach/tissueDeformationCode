#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 12:17:11 2022

@author: Javiera Jilberto Vallejos
"""

from PIL import Image
from skimage.transform import warp
from skimage.color import rgb2gray
import numpy as np
from numpy import asarray
import matplotlib.pyplot as plt

# Load image and converted to grayscale



def transformImage(imageLocation):
    image = Image.open(imageLocation)
    image = rgb2gray(np.array(image)[:,:,0:3])

    def xyStretch(originalImage):
        # Define transform
        T = np.array([[1.2, 0.1],
                      [0.1, 1]])

        # Compute original image center and deform image center
        im_center = np.asarray(image.T.shape)/2
        A = T@im_center - im_center

        # Define inverse of the transform
        W = np.eye(3)
        W[0:2,0:2] = np.linalg.inv(T)
        W[0:2,2] = A

        # Warp image
        wim = warp(image, W)
        return wim
    
    def imageToArray(imageLoc):
        image = Image.open(imageLoc+'.png')
        image = rgb2gray(np.array(image)[:,:,0:3])
        arrayImg = asarray(image)
        print(arrayImg)
        return arrayImg

    warpedImage = xyStretch(image)
    
    np.save('originalArray.npy', image)
    np.save('transformedArray.npy', warpedImage)
      
    # Visualize
    plt.figure(1, clear=True, frameon=False)
    plt.imshow(image, vmin=0, vmax=1)
    plt.axis('off')
    plt.savefig('initialImage.png')

    plt.figure(2, clear=True, frameon=False)
    plt.imshow(warpedImage, vmin=0, vmax=1)
    plt.axis('off')
    plt.savefig('transformedImage.png')

    
transformImage('newTestImage.png')



