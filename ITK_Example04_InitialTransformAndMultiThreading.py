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
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Import Images
fixed_image = itk.imread('transformingTestImage/newTestImage.png', itk.F)
moving_image = itk.imread('transformingTestImage/newTransformedImage.png', itk.F)
# fixed_image = itk.imread('data/align0.1LAP_05CROPPED_30.png', itk.F)
# moving_image = itk.imread('data/align0.1LAP_05CROPPED_55.png', itk.F)

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
plt.axis('off')

plt.imshow(result_image)

itk.imwrite(result_image, 'output/result_image_11_7.nii')

#########Deformation Field#########
deformation_field = itk.transformix_deformation_field(moving_image, result_transform_parameters)

array = itk.GetArrayFromImage(fixed_image)
file = open('output/outputArrayOG.txt', 'w')
file.write(" ".join(str(x) for x in array))
file.close()

deformation_field2 = np.asarray(deformation_field).astype(np.float32)

# # Write the raw data to a file
# # print(deformation_field)
# # file = open('output/outputDeformationField.txt', 'w')
# # file.write(" ".join(str(x) for x in deformation_field))
# # file.close()

# #Plot images
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


# plt.colorbar()
