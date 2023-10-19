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
from skimage.transform import PiecewiseAffineTransform, warp



def createImage(image_fldr, output_fldr):
    image = imread(image_fldr, as_gray=True)
    
    a = 1.05

    def stretch(image):

        i = np.arange(image.shape[1])
        j = np.arange(image.shape[0])
        I, J = np.meshgrid(i,j)
    
        im_size = J.shape
        im_center = np.asarray(im_size)/2
        
        # Mapping function
        I_center = I-im_center[0]
        J_center = J-im_center[1]
        
        # Transformed array
        T = np.zeros([2, *I_center.shape])
        
        
        IJ = np.zeros([2, *I.flatten().shape])
        T2 = np.zeros([2, *I.flatten().shape])
        IJ[0] = I.flatten()
        IJ[1] = J.flatten()

        I_disp = a*I_center
        J_disp = J_center*0
        
        # print(I_disp.shape)
        
        T[0] = I_disp + (I_center + im_center[0])
        T[1] = J_disp + (J_center + im_center[1])
        
        T2[0] = T[0].flatten()
        T2[1] = T[1].flatten()
        
        tform = PiecewiseAffineTransform()
        tform.estimate(IJ.T, T2.T)

        #create analytical deformation gradient
        analytical_gradient = []
        for i in range (T[0].shape[0]):
            row = []
            for j in range (T[0].shape[1]):
                sub_matrix = [[3, 0], [0, 1]]
                row.append(np.asarray(sub_matrix))
            analytical_gradient.append(np.asarray(row))
        
        out = warp(image, tform.inverse)
        return out, analytical_gradient, I_disp, J_disp, T
        
    warpedImage, analytical_gradient, I_disp, J_disp, T  = stretch(image)
    imsave(output_fldr, warpedImage)
    
    

    # Visualize
    fig, axs = plt.subplots(2, 1, num=1, clear=True)
    axs[0].imshow(image)
    axs[1].imshow(warpedImage)
    
    #analytical displacement
    return analytical_gradient, I_disp, J_disp, T
    



