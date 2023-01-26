
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 16:43:42 2023

@author: Javiera Jilberto Vallejos
"""

from skimage.transform import warp
from skimage.io import imread
from skimage.color import rgb2gray
import numpy as np
import matplotlib.pyplot as plt

# Load image and converted to grayscale
image = imread('newTestImage.png', as_gray=True)
# image = rgb2gray(np.array(image)[:,:,0:3])


i = np.arange(image.shape[1])
j = np.arange(image.shape[0])
I, J = np.meshgrid(i,j)

im_size = J.shape
im_center = np.asarray(im_size)/2

# Mapping function
I_center = J-im_center[0]
J_center = I-im_center[1]

# Transformed array
T = np.zeros([2, *I_center.shape])

# I_disp = J_center**2*1e-3 - 10
# J_disp = I_center*1e-1


I_disp = np.cos(J_center/20)*10
J_disp = J_center/1e3

T[0] = I_disp + (I_center + im_center[0])
T[1] = J_disp + (J_center + im_center[1])

warp_image = warp(image, T)


fig, axs = plt.subplots(2, 2, num=1, clear=True)
axs[0,0].imshow(image)
axs[0,1].imshow(warp_image)
axs[1,0].imshow(I_disp)
axs[1,1].imshow(J_disp)

