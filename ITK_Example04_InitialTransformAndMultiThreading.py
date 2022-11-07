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

# In[1]:


import itk
import matplotlib.pyplot as plt

# In[7]:


# Import Images
fixed_image = itk.imread('data/newTestImage.png', itk.F)
moving_image = itk.imread('data/newTransformedImage.png', itk.F)

# Import Default Parameter Map
parameter_object = itk.ParameterObject.New()
resolutions = 3
parameter_map_rigid = parameter_object.GetDefaultParameterMap('rigid',3)
parameter_object.AddParameterMap(parameter_map_rigid)

# For the bspline default parameter map, an extra argument can be specified that define the final bspline grid spacing in physical space. 
parameter_map_bspline = parameter_object.GetDefaultParameterMap("bspline", resolutions, 20.0)
parameter_object.AddParameterMap(parameter_map_bspline)


# .. and/or load custom parameter maps from .txt file
parameter_object.AddParameterFile('data/parameters_BSpline.txt')

# Call registration function with initial transfrom and number of threads
result_image, result_transform_parameters = itk.elastix_registration_method(
    fixed_image, moving_image,
    parameter_object=parameter_object,
    log_to_console=False)

# Update filter object (required)
# elastix_object.UpdateLargestPossibleRegion()

# # Results of Registration
# result_image = elastix_object.GetOutput()
# result_transform_parameters = elastix_object.GetTransformParameterObject()
array = itk.GetArrayFromImage(result_image)
plt.imshow(result_image)
plt.savefig(r'transform11_2.png')
print(array)


# Pointing from a transform parameter file to the path of a second initial transform parameter file is supported from the 0.7.0 release of ITKElastix.
