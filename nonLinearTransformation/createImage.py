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
from strainOperations import display_strain


def createImage(image_fldr, output_fldr):
    image = imread(image_fldr, as_gray=True)

    a = 20
    b = 10

    def stretch(image):
        i = np.arange(image.shape[0])
        j = np.arange(image.shape[1])
        I, J = np.meshgrid(i,j)
        I = I.T
        J = J.T

        im_size = J.shape

        im_center = np.asarray(im_size)/2

        # Mapping function
        I_center = I-im_center[0]
        J_center = J-im_center[1]

        # Transformed array
        T = np.zeros([2, *I_center.shape])

        I_disp =(np.cos(J_center/a)*b)
        J_disp = np.zeros_like(I_disp)

        T[0] = I_center - (b*np.cos(J_center/a)) + im_center[0]
        T[1] = J_center + im_center[1]

        #create analytical deformation gradient
        analytical_gradient = []
        for i in range (T[0].shape[0]):
            row = []
            for j in range (T[0].shape[1]):
                alpha = (b/a)*np.sin((j-im_center[1])/a)
                sub_matrix = [[1, alpha], [0, 1]]

                row.append(np.asarray(sub_matrix))
            analytical_gradient.append(np.asarray(row))

        warp_image = warp(image, T)

        E = np.zeros_like(analytical_gradient)
        alpha = (b/a)*np.sin(J_center/a)
        E[:,:,0,1] = ((b/a)*np.sin(J_center/a))/2
        E[:,:,1,0] = ((b/a)*np.sin(J_center/a))/2
        E[:,:,1,1] = (((b/a)*np.sin(J_center/a))**2)/2


        return warp_image, analytical_gradient, I_disp, J_disp, E

    warpedImage, analytical_gradient, I_disp, J_disp, E  = stretch(image)
    imsave(output_fldr, warpedImage)

    # Visualize
    fig, axs = plt.subplots(2, 1, num=1, clear=True)
    axs[0].imshow(image)
    axs[1].imshow(warpedImage)
    # axs[0].colorbar()
    #analytical displacement
    return analytical_gradient, I_disp, J_disp, E




