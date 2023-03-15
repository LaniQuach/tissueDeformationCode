# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:26:00 2023

@author: laniq
"""

from skimage import io
import os
import matplotlib.pyplot as plt

# Load .tif stack
def getFullStack(tissueName):
    
    im = io.imread(tissueName + '.tif')
    # To save all frames (it will create a bunch of files, so might want to create a specific folder for this)
    ext = '.png'
    for frame in range(0, im.shape[0]):
        completeName = tissueName + '_%i' % frame + ext
        #completeName = os.path.join(path, filename)
        io.imsave(completeName, im[frame])
    
def getIndividualImages(tissueName, frame):
    # To save a frame
    im = io.imread(tissueName + '.tif')

    ext = '.png'
    completeName = tissueName + '_%i' % frame + ext
    io.imsave(completeName, im[frame])



