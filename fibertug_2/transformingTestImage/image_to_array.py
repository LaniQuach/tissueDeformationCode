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


image = imread('A_0.41_01-14-2_rc_ds_0.png', as_gray=True)
image_2 = imread('A_0.41_01-14-2_rc_ds_25.png', as_gray=True)

print(np.max(np.asarray(image)))
print(np.max(np.asarray(image_2)))

image_array = np.asarray(image).astype(np.float32)
image_2_array = np.asarray(image_2).astype(np.float32)

mask_1 = np.asarray(imread('Mask_0.png', as_gray=True))
mask_2 = np.asarray(imread('Mask_25.png', as_gray=True))

image_array[mask_1==0] = 0
image_2_array[mask_2==0] = 0

np.save('originalArray.npy', image_array)
np.save('transformedArray.npy', image_2_array)
  
# Visualize
fig, axs = plt.subplots(2, 1, num=1, clear=True)
axs[0].imshow(image_array)
axs[0].axis('off')
axs[0].set_title('Relaxed')

axs[1].imshow(image_2_array)
axs[1].axis('off')
axs[1].set_title('Contracted')

plt.savefig('transformedImage.png')

# transformImage('align0.1LAP_05CROPPED_31.png', 'align0.1LAP_05CROPPED_45.png')



