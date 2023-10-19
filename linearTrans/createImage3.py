# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:15:51 2023

@author: laniq
"""

from skimage.transform import warp
from skimage.io import imread
from skimage.io import imsave
import numpy as np
import matplotlib.pyplot as plt

def createImage(image_fldr, output_fldr):
    image = imread(image_fldr, as_gray=True)
    
    a = 1

    def stretch(image):
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
        
        I_disp = 1/(I_center*a)
        J_disp = J_center*0
        
        # print(I_disp.shape)
        
        T[0] = I_disp + (I_center + im_center[0]))
        T[1] = J_disp + (J_center + im_center[1])

        #create analytical deformation gradient
        analytical_gradient = []
        for i in range (T[0].shape[0]):
            row = []
            for j in range (T[0].shape[1]):
                sub_matrix = [[a+1,0], [0,1]]
                row.append(np.asarray(sub_matrix))
            analytical_gradient.append(np.asarray(row))

        warp_image = warp(image, T)
        
        return warp_image, analytical_gradient, I_disp, J_disp, T
        
    warpedImage, analytical_gradient, I_disp, J_disp, T  = stretch(image)
    imsave(output_fldr, warpedImage)
    
    

    # Visualize
    fig, axs = plt.subplots(2, 1, num=1, clear=True)
    axs[0].imshow(image)
    axs[1].imshow(warpedImage)
    
    #analytical displacement
    return analytical_gradient, I_disp, J_disp, T
    



