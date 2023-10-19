# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:57:07 2023

@author: laniq
"""



from skimage.transform import warp
from skimage.io import imread
from skimage.io import imsave
import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import PiecewiseAffineTransform, warp
from PIL import Image

def createImage(image_fldr, output_fldr):
    # image = imread(image_fldr, as_gray=True)
    image = Image.open(image_fldr).convert('L')
    a = 2

    def stretch(image):
        width, height = image.size
        resized = image.resize((width*a, height))
        width2, height2 = resized.size
 
        #crop transformed image
        left = (width2-width)/2
        right = left+width
        top = 0
        bottom = height
        resizedCrop = resized.crop((left, top, right, bottom))
        
        #create analytical deformation gradient
        analytical_gradient = []
        for i in range (np.asarray(image).shape[0]):
            row = []
            for j in range (np.asarray(image).shape[1]):
                sub_matrix = [[3, 0], [0, 1]]
                row.append(np.asarray(sub_matrix))
            analytical_gradient.append(np.asarray(row))
        
        
        return resizedCrop, analytical_gradient
        
    warpedImage, analytical_gradient  = stretch(image)
    imsave(output_fldr, np.asarray(warpedImage))
    
    

    # Visualize
    fig, axs = plt.subplots(2, 1, num=1, clear=True)
    axs[0].imshow(image)
    axs[1].imshow(warpedImage)
    
    #analytical displacement
    return analytical_gradient
    



