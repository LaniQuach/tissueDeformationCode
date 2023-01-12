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
import matplotlib.pyplot as plt

# Load image and converted to grayscale
image = Image.open('newTestImage.png')
image = rgb2gray(np.array(image)[:,:,0:3])

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

# Visualize
plt.figure(1, clear=True)
plt.imshow(image)

plt.figure(2, clear=True)
plt.imshow(wim)
