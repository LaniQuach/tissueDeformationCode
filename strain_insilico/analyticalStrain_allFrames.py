# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 09:27:43 2023

@author: laniq
"""
import os
import numpy as np
from matplotlib import pyplot as plt
import cheartio as chio
from skimage.io import imread
from matplotlib.tri import Triangulation, LinearTriInterpolator
from strainOperations import secondDeriv, calculate_defGradient, calculate_strain, display_strain, display_strainC, plot_deriv, calculate_strain_withC
from tqdm import tqdm


sim_fldr = 'sim_data/'
img_fldr = 'warped_image/images/'
pix_vid = 0.908
num_frames = 51
fldr = 'strain_calc/'
out_fldr = fldr + 'out/'
mask1 = fldr + 'Mask_1.png'

if not os.path.exists(out_fldr): os.mkdir(out_fldr)

parameterFileName = fldr + 'parameters_BSpline.txt'
original_image = imread(img_fldr + 'warped_0' + '.png', as_gray = True)

original_array = np.asarray(original_image).astype(np.float32).T
mask_1 = np.asarray(imread(mask1, as_gray=True))
mask_1 = mask_1[9:-10]
Mask_1 = np.zeros((244, 420))
Mask_1[0:mask_1.shape[0],11:mask_1.shape[1]+11] = mask_1
# mask_1 = np.resize(mask_1, (244, 420))
plt.imshow(Mask_1)


i = np.arange(original_image.shape[0])
j = np.arange(original_image.shape[1])
I,J = np.meshgrid(i,j)

ijk = np.vstack([I.flatten(), J.flatten()]).T   # pixel coordinates

xyz_quad, ien_quad, _ = chio.read_mesh(sim_fldr + 'mesh/tissue_quad')
xyz_quad = xyz_quad*1000/pix_vid          # Transform to pixel space
tri_quad = Triangulation(xyz_quad[:,0], xyz_quad[:,1], triangles=ien_quad[:,0:3])

disp = chio.read_dfile(sim_fldr + 'out/' + 'U-%i' % 25 + '.D')
disp = disp*1000/pix_vid

image = imread(img_fldr + 'warped_%i' % 25 + '.png')
image_array = np.asarray(image).astype(np.float32).T

xyz, ien, _ = chio.read_mesh(sim_fldr + 'mesh/tissue')
xyz = xyz*1000/pix_vid          # Transform to pixel space

# Translate mesh
tx = 18/pix_vid
ty = 12/pix_vid
xyz[:,0] += tx
xyz[:,1] += ty

trimesh = Triangulation(xyz[:,0],xyz[:,1],triangles=ien)

v_min1 = -0.35
v_max1 = 0.7

v_min2 = -0.3
v_max2 = 0.4

v_min4 = -0.3
v_max4 = 0.4

# interp_func2 = LinearTriInterpolator(trimesh, strains_analyt[:,1])
# img_frho2 = interp_func2(ijk[:,1], ijk[:,0]).data.reshape(I.shape).T
# img_frho2 = np.flipud(img_frho2)
# # img_frho2 = img_frho2[2:-2, 2:-2]

# interp_func3 = LinearTriInterpolator(trimesh, strains_analyt[:,3])
# img_frho3 = interp_func3(ijk[:,1], ijk[:,0]).data.reshape(I.shape).T
# img_frho3 = np.flipud(img_frho3)
# img_frho3 = img_frho3[2:-2, 2:-2]
sd_boxList = np.load("subdomainBoxes.npy")
coordinates = np.zeros((len(sd_boxList), 3))
strainValues = np.zeros((len(sd_boxList), num_frames))

# strain in ecc
for ts in tqdm(range(1,num_frames)):
    #strains using analytical exact
    C = chio.read_dfile(sim_fldr + 'out/RCauchyGreen-%i' % (ts+1) + '.D')
    
    #use a for loop and go over the C array
    # each line has 4 values, turn that into a 2x2 matrix
    strains_analyt = np.zeros(C.shape)
    
    C_array = np.asarray(C)
    for i in range (C.shape[0]):
        x = np.asarray(np.reshape(C[i], (2,2)))
        strains_analyt[i] = calculate_strain_withC(x)

    interp_func1 = LinearTriInterpolator(trimesh, strains_analyt[:,0])
    img_frho1 = interp_func1(ijk[:,1], ijk[:,0]).data.reshape(I.shape).T
    # img_frho1 = np.flipud(img_frho1)

    # plt.figure()
    # plt.imshow(img_frho1)
    # plt.scatter([upper_left_y, upper_right_y, lower_left_y, lower_right_y], [upper_left_x, upper_right_x, lower_left_x, lower_right_x])
    # plt.scatter(coordinates[sd][0], coordinates[sd][1])
    plt.show()
    # calculate the average over all points in the sd
    for sd in range(0, len(sd_boxList)):
        # Define the coordinates of the subdomain box
        upper_left_x = sd_boxList[sd][0][0]
        upper_left_y = sd_boxList[sd][0][1]
        upper_right_x = sd_boxList[sd][1][0]
        upper_right_y = sd_boxList[sd][1][1]
        lower_left_x = sd_boxList[sd][3][0]
        lower_left_y = sd_boxList[sd][3][1]
        lower_right_x = sd_boxList[sd][2][0]
        lower_right_y = sd_boxList[sd][2][1]
        
        coordinates[sd][1] = (upper_left_x + upper_right_x + lower_left_x + lower_right_x) / 4.0
        coordinates[sd][0] = (upper_left_y + upper_right_y + lower_left_y + lower_right_y) / 4.0
        
        
        # Initialize variables to store the sum of values and count of valid pixels
        sum_of_values = 0
        count_of_valid_pixels = 0
        
        # Loop through the region defined by the coordinates
        for y in range(upper_left_x, lower_left_x + 1):
            for x in range(upper_left_y, upper_right_y + 1):
                
                # Check if the pixel coordinates are within the region
                if not np.isnan(img_frho1[y, x]):
                    # Get the pixel value from img_frho1
                    pixel_value = img_frho1[y, x]
                    # Add the pixel value to the sum
                    sum_of_values += pixel_value
                    # Increment the count of valid pixels
                    count_of_valid_pixels += 1
        
        # Calculate the average
        if count_of_valid_pixels > 0:
            average_value = sum_of_values / count_of_valid_pixels
        else:
            average_value = 0
            
        strainValues[sd][ts] = average_value
        coordinates[sd][2] = average_value
        

        
    # Now, 'average_value' contains the average of all valid non-NaN values within the subdomain box
    image = imread(img_fldr + 'warped_%i' % ts + '.png', as_gray = True)

    plt.figure()
    # plt.imshow(image)
    plt.imshow(image, cmap='Greys')
    plt.scatter(coordinates[:,0], coordinates[:,1], c = coordinates[:,2], s = 50, cmap = plt.cm.RdBu, vmin=-0.09, vmax=0.0)
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.title("Analytical Strain frame: %i" %ts)
    plt.show()

np.save("Subdomain_ECC_analytical_finite.npy", strainValues)

opticalFlow = np.load("Subdomain_ECC_analytical.npy")


