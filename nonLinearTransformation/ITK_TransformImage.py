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
from numpy import asarray

from mpl_toolkits.axes_grid1 import make_axes_locatable

#load both arrays
fixed_array = np.load('transformingTestImage/originalArray.npy')
moving_array = np.load('transformingTestImage/transformedArray.npy')

fixed_array = fixed_array.astype(np.float32)
moving_array = moving_array.astype(np.float32)

fixed_image = itk.GetImageFromArray(fixed_array)
moving_image = itk.GetImageFromArray(moving_array)

# show images
plt.figure(frameon=False)
plt.imshow(fixed_image, vmin = 0, vmax = 1)
plt.axis('off')
plt.savefig('output/originalImage.png')

plt.figure(frameon=False)
plt.imshow(moving_image, vmin = 0, vmax = 1)
plt.axis('off')
plt.savefig('output/stretchedImage.png')

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
result_array = itk.GetArrayFromImage(result_image)

result_transform_parameters = elastix_object.GetTransformParameterObject()
np.save('output/elastixResultImage_Array.npy', result_array)

plt.figure(frameon=False)
plt.imshow(result_image, vmin = 0, vmax = 1, aspect = 'auto')
plt.axis('off')
plt.savefig('output/elastixResults.png')
# itk.imwrite(result_image, 'output/result_image.nii')

######### Deformation Field #########
deformation_field = itk.transformix_deformation_field(moving_image, result_transform_parameters)
defArray = itk.GetArrayFromImage(deformation_field).astype(float)*0.908

np.save('output/analytical_x.npy', defArray[:,:,0])
np.save('output/analytical_y.npy', defArray[:,:,1])

fig, axs = plt.subplots(1, 2, sharey=True, figsize=[30,30])
im3 = axs[1].imshow(defArray[:,:,0])
#, vmin = -15, vmax = 15
divider = make_axes_locatable(axs[1])
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar = fig.colorbar(im3, cax=cax, orientation='vertical');
cbar.set_label('displacement (pixels)', fontsize = 25)
cbar.ax.tick_params(labelsize=30)


im2 = axs[0].imshow(defArray[:,:,1]*-1)
#, vmin = -5, vmax = 5
divider1 = make_axes_locatable(axs[0])
cax1 = divider1.new_vertical(size='5%', pad=0.2, pack_start = True)
fig.add_axes(cax1)
cbar1 = fig.colorbar(im2, cax = cax1, orientation = 'horizontal')
cbar1.set_label('displacement (pixels)', fontsize = 25)
cbar1.ax.tick_params(labelsize=25)


axs[0].axis('off')
axs[1].axis('off')
axs[0].set_title('Displacement Field Y', fontsize=30)
axs[1].set_title('Displacement Field X', fontsize=30)


plt.savefig('output/' + 'displacement' + '.png', dpi = 200)


#as a note next time save the above file with the colorbar being the same for both

