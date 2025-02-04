# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 13:36:39 2023

@author: laniq
"""
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from ITK_TransformImage import *
from tqdm import tqdm

data = np.load('data/post_displacement.npy', allow_pickle=True).item()
data_emma = data['emma'][0:51]
data_sam = data['sam'][0:51]

beat = 1
pos = {}
pos['x'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_col.txt')
pos['y'] = np.loadtxt('comparisionPlots' + '/beat%i' % beat + '_row.txt')

point = 36
emmas_x = np.round(pos['x'][point,:],0).astype(int)
emmas_y = np.round(pos['y'][point,:],0).astype(int)

query_point = np.array([137, 186])

original_image = imread('data/fibers_2' + '/A_0.41_01-14-2_rc_ds_%i' % 0 + '.png', as_gray=True)
plt.figure()
plt.plot(emmas_x[0], emmas_y[0], marker='.', color="white")
plt.imshow(original_image, vmin = 0, vmax = 10280)
plt.axis('off')
plt.show()
plt.savefig('comparisionPlots/locationOfPost.png')


original_array = np.asarray(original_image).astype(np.float32).T
elastix_x = np.empty([51])
elastix_x[0] = 0

for ts in tqdm(range(1,51)):
    image = imread('data/fibers_2' + '/A_0.41_01-14-2_rc_ds_%i' % ts + '.png', as_gray=True)
    image_array = np.asarray(image).astype(np.float32).T

    resultTransformParameters = elastix_transformation(original_array, image_array)
    disp_array = displacement_field_elastix(resultTransformParameters, image_array)

    elastix_x[ts] = -np.abs(disp_array[query_point[1], query_point[0], 0])


x = range(0,51)
plt.figure()
plt.plot(x, data_emma, color='blue', label = 'Optical Flow')
plt.plot(x, elastix_x*-1, color='magenta', label = 'Elastix')
plt.plot(x, data_sam, color='green', label = 'sam')
plt.title('Post Displacement')
# axs[1].annotate('R2=%2.3f' % r2_y, (0.05,0.9), xycoords='axes fraction')
plt.xlabel("frame")
plt.ylabel("displacement")
plt.legend()
plt.show()
plt.savefig('comparisionPlots/PostDisplacementPlot.png')




