#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 16:43:42 2023

@author: Javiera Jilberto Vallejos
"""

from skimage.transform import warp
from skimage.io import imread
from skimage.io import imsave
import numpy as np
import matplotlib.pyplot as plt


def createImage(image_fldr, output_fldr):
    image = imread(image_fldr, as_gray=True)
    
    a = 20
    b = 10
    c = 1e3

    def stretch(image):
        # warp_image = np.zeros(image.shape, dtype = image.dtype)
        
        # rows, cols = image.shape
        # for i in range(rows):
        #     for j in range(cols):
        #         offset_x = int(j/c)
        #         offset_y = int(np.cos(j/a)*b)
        #         if i+offset_y < rows:
        #             warp_image[i,j] = image[(i + offset_y), (j + offset_x)]
        #         else:
        #             warp_image[i,j] = 0
                    
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
        
        I_disp = np.cos(J_center/a)*b
        J_disp = J_center/c
        
        T[0] = I_disp + (I_center + im_center[0])
        T[1] = J_disp + (J_center + im_center[1])

        #create analytical deformation gradient
        analytical_gradient = []
        for i in range (T[0].shape[0]):
            row = []
            for j in range (T[0].shape[1]):
                sub_matrix = [[1, (1/c)], [0, ((b/a)*np.sin(j/a))+1]]
                row.append(np.asarray(sub_matrix))
            analytical_gradient.append(np.asarray(row))

        warp_image = warp(image, T)
        return warp_image, analytical_gradient
        
    warpedImage, analytical_gradient = stretch(image)
    imsave(output_fldr, warpedImage)

    # Visualize
    fig, axs = plt.subplots(2, 1, num=1, clear=True)
    axs[0].imshow(image)
    axs[1].imshow(warpedImage)
    
    #analytical displacement
    return analytical_gradient
    



