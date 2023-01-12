#!/usr/bin/env python
# coding: utf-8

# ## 4. Image Registration with initial transform and/or multiple threads

# In this notebook 2 other options of the elastix algorithm are shown: initial transformation and multithreading.
# They're shown together just to reduce the number of example notebooks and 
# thus can be used independently as well as in combination with whichever other functionality
# of the elastix algorithm. 
# 
# Initial transforms are transformations that are done on the moving image before the registration is started.
# 
# Multithreading spreaks for itself and can be used in similar fashion in the transformix algorithm.
# 
# 

# ### Registration

import itk
import SimpleITK as sitk

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Import Images
# fixed_image = itk.imread('transformingTestImage/newTestImage.png', itk.F)
# moving_image = itk.imread('transformingTestImage/newTransformedImage.png', itk.F)

fixed_array = np.load('transformingTestImage/original.npy')
moving_array = np.load('transformingTestImage/transformed.npy')

fixed_array = fixed_array.astype(np.float32)
moving_array = moving_array.astype(np.float32)

fixed_image = itk.GetImageFromArray(fixed_array)
moving_image = itk.GetImageFromArray(moving_array)

plt.figure(frameon=False)
plt.imshow(fixed_array, vmin=0, vmax=255, aspect='auto')
plt.axis('off')
plt.savefig('output/ogImage1.png')

plt.figure(frameon=False)
plt.imshow(moving_array, vmin=0, vmax=255, aspect='auto')
plt.axis('off')
plt.savefig('output/movingImage.png')

# Import Default Parameter Map
parameter_object = itk.ParameterObject.New()

# Load custom parameter maps from .txt file
parameter_object.AddParameterFile('data/parameters_BSpline.txt')

# Load Elastix Image Filter Object
elastix_object = itk.ElastixRegistrationMethod.New(fixed_image, moving_image)
elastix_object.SetParameterObject(parameter_object)

# Set additional options
elastix_object.SetLogToConsole(False)

# Update filter object (required)
elastix_object.UpdateLargestPossibleRegion()

# Results of Registration
result_image = elastix_object.GetOutput()
result_transform_parameters = elastix_object.GetTransformParameterObject()

plt.figure(frameon=False)
plt.imshow(result_image, vmin=0, vmax=255, aspect='auto')
plt.axis('off')
plt.savefig('output/elastixResults_1_12.png')

# plt.figure()
# plt.imshow(fixed_image)
# plt.axis('off')

# plt.savefig('output/ogImage.png')

itk.imwrite(result_image, 'output/result_image_11_7.nii')

#########Deformation Field#########
deformation_field = itk.transformix_deformation_field(moving_image, result_transform_parameters)

array = itk.GetArrayFromImage(fixed_image)
file = open('output/outputArrayOG.txt', 'w')
file.write(" ".join(str(x) for x in array))
file.close()

deformation_field2 = np.asarray(deformation_field).astype(np.float32)

# Write the raw data to a file
# print(deformation_field)
# file = open('output/outputDeformationField.txt', 'w')
# file.write(" ".join(str(x) for x in deformation_field))
# file.close()

#Plot images
fig, axs = plt.subplots(1, 2, sharey=True, figsize=[30,30])
plt.figsize=[100,100]
im2 = axs[0].imshow(deformation_field[:,:,1])
divider = make_axes_locatable(axs[0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2, cax=cax, orientation='vertical');
axs[0].set_title('Deformation Field X', fontsize=30)

im3 = axs[1].imshow(deformation_field[:,:,0])
axs[1].set_title('Deformation Field Y', fontsize=30)
divider = make_axes_locatable(axs[1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im3, cax=cax, orientation='vertical');
plt.savefig('output/deformationField.png')


# plt.colorbar()
