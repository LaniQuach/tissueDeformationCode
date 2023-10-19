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

def imageToArray(imageFileName1, imageFileName2, mask1=None, mask2=None):
    """
    image_to_array transform png images to arrays to be used in elastix transformations

    :param imageFileName1: original image of the image before contraction/movement
    :param imageFileName2: image during contraction/movement
    :param mask1, mask2: masks for image1 and image2
    :return: the array versions of image1 and image 2 with the masks applied as a list of two values
             list[0] = imageArray1 and list[1] = imageArray2
    """
    image = imread(imageFileName1, as_gray=True)
    image_2 = imread(imageFileName2, as_gray=True)

    #transforms images into arrays
    image_array = np.asarray(image).astype(np.float32)
    image_2_array = np.asarray(image_2).astype(np.float32)

    #converts masks into arrays
    if mask1 is not None:
        mask_1 = np.asarray(imread(mask1, as_gray=True))
        image_array[mask_1==0] = 0
    if mask2 is not None:
        mask_2 = np.asarray(imread(mask2, as_gray=True))
        image_2_array[mask_2==0] = 0

    return [image_array, image_2_array]

def display_bothImages(image_array, image_2_array, name, save):
    # Visualize
    fig, axs = plt.subplots(2, 1, num=1, clear=True)
    axs[0].imshow(image_array, cmap = 'Greys')

    axs[0].axis('off')
    axs[0].set_title('Target')

    axs[1].imshow(image_2_array, cmap = 'Greys')
    axs[1].axis('off')
    axs[1].set_title('Stretched')

    if save:
        plt.savefig(name + '.png', dpi = 210)



