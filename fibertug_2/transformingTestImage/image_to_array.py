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

np.save('originalArray.npy', image)
np.save('transformedArray.npy', image_2)
  
# Visualize
fig, axs = plt.subplots(2, 1, num=1, clear=True)
axs[0].imshow(image)
axs[0].axis('off')
axs[0].set_title('Relaxed')

axs[1].imshow(image_2)
axs[1].axis('off')
axs[1].set_title('Contracted')

plt.savefig('transformedImage.png')

# transformImage('align0.1LAP_05CROPPED_31.png', 'align0.1LAP_05CROPPED_45.png')



