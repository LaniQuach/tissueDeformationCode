# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 23:14:28 2023

@author: laniq
"""

from ITK_TransformImage import displacement_field_elastix, display_save_displacement
from image_to_array import imageToArray, display_bothImages
from compareImages import compareStrain
from strainOperations import secondDeriv, calculate_defGradient, calculate_strain, display_strain, display_strainC, plot_deriv, calculate_strain_withC
import os
import numpy as np
from matplotlib import pyplot as plt
import cheartio as chio
from skimage.io import imread
import math
from matplotlib.tri import Triangulation, LinearTriInterpolator

sim_fldr = 'sim_data/'
img_fldr = 'warped_image/images/'
pix_vid = 0.908
fldr = 'strain_calc/'
out_fldr = fldr + 'out/'
mask1 = fldr + 'tissue_mask.png'
if not os.path.exists(out_fldr): os.mkdir(out_fldr)

parameterFileName = fldr + 'parameters_BSpline.txt'

# Read image
ts = 25
targetImageName = img_fldr + 'warped_0' + '.png'
stretchedImageName = img_fldr + 'warped_%02i' % ts + '.png'

warp_img = imread(img_fldr + 'warped_%02i' % ts + '.png')

imageArrays = imageToArray(targetImageName, stretchedImageName)

display_bothImages(imageArrays[0], imageArrays[1], out_fldr + 'originalImages', True)

#Elastix functionality
displacementField = displacement_field_elastix(imageArrays[0], imageArrays[1], parameterFileName, True, out_fldr, pix_vid, mask1) #u
# display_save_displacement(displacementField, tri_quad, disp, out_fldr + 'dispField', True)

trackingpoints = np.load('tracking_points.npy')
xcoords = trackingpoints[:,0,0]
xcoords = [ int(x) for x in xcoords ]
ycoords = trackingpoints[:,0,1]
ycoords = [ int(x) for x in ycoords ]

numpoints = ycoords.shape[0]
elast_tracked0 = np.zeros((numpoints,2)) #x dir
elast_tracked1 = np.zeros((numpoints,2)) #y dir

for xx in range (0, numpoints):
    #if coordinate point at xx is in the displacement field in both the x and y direction
    if(xcoords[xx] < displacementField.shape[1]) and (ycoords[xx] < displacementField.shape[0]): 
        elast_tracked0[xx,0] = ycoords[xx]
        elast_tracked1[xx,0] = xcoords[xx]
        
        elast_tracked0[xx,1] = ycoords[xx] + displacementField[ycoords[xx], xcoords[xx], 1]  
        elast_tracked1[xx,1] = xcoords[xx] + displacementField[ycoords[xx], xcoords[xx], 0]
        
        # print("xdir", elast_tracked0[xx,1])
        # print("ydir", elast_tracked1[xx,1])
        
    
# print(analytic_tracked0[20])
np.savetxt('xdir_elastixDisp.txt', elast_tracked0)
np.savetxt('ydir_elastixDisp.txt', elast_tracked1)

# # #strains using elastix
# # derivatives_elastix = secondDeriv(displacementField[:,:,1], displacementField[:,:,0])
# # # plot_deriv(derivatives_elastix, out_fldr, False)
# # defGradient_elastix = calculate_defGradient(derivatives_elastix) #gradient u
# # strains_elastix = calculate_strain(defGradient_elastix)
# # # display_strain(strains_elastix, out_fldr + 'elastix_', True)

# #strains using analytical exact
# C = chio.read_dfile(sim_fldr + 'out/U-%i' % ts + '.D')
# #use a for loop and go over the C array
# # each line has 4 values, turn that into a 2x2 matrix
# # strains_analyt = np.zeros(C.shape)

# C_array = np.asarray(C)

# for i in range (C.shape[0]):
#     x = np.asarray(np.reshape(C[i], (2,2)))
#    strains_analyt[i] = calculate_strain_withC(x)


# xyz, ien, _ = chio.read_mesh(sim_fldr + 'mesh/tissue')
# xyz = xyz*1000/pix_vid          # Transform to pixel space

# # Translate mesh
# tx = 18/pix_vid
# ty = 12/pix_vid
# xyz[:,0] += tx
# xyz[:,1] += ty

# trimesh = Triangulation(xyz[:,0],xyz[:,1],triangles=ien)

# v_min1 = -0.35
# v_max1 = 0.7

# v_min2 = -0.3
# v_max2 = 0.4

# v_min4 = -0.3
# v_max4 = 0.4

# display_strainC(trimesh, strains_analyt[:,0], strains_elastix[1:-1,1:-1,0,0], out_fldr + 'analytical_', v_min1, v_max1, 'xx', True)
# display_strainC(trimesh, strains_analyt[:,1], strains_elastix[1:-1,1:-1,0,1], out_fldr + 'analytical_', v_min2, v_max2, 'xy' , True)
# display_strainC(trimesh, strains_analyt[:,3], strains_elastix[1:-1,1:-1,1,1], out_fldr + 'analytical_', v_min4, v_max4, 'yy', True)

# interp_func1 = LinearTriInterpolator(trimesh, strains_analyt[:,0])
# img_frho1 = interp_func1(ijk[:,1], ijk[:,0]).data.reshape(I.shape).T
# img_frho1 = np.flipud(img_frho1)
# img_frho1 = img_frho1[2:-2, 2:-2]

# interp_func2 = LinearTriInterpolator(trimesh, strains_analyt[:,1])
# img_frho2 = interp_func2(ijk[:,1], ijk[:,0]).data.reshape(I.shape).T
# img_frho2 = np.flipud(img_frho2)
# img_frho2 = img_frho2[2:-2, 2:-2]

# interp_func3 = LinearTriInterpolator(trimesh, strains_analyt[:,3])
# img_frho3 = interp_func3(ijk[:,1], ijk[:,0]).data.reshape(I.shape).T
# img_frho3 = np.flipud(img_frho3)
# img_frho3 = img_frho3[2:-2, 2:-2]

# compareStrain(strains_elastix[1:-1,1:-1,0,0], img_frho1, v_min1, v_max1, "xx", out_fldr)
# compareStrain(strains_elastix[1:-1,1:-1,0,1], img_frho2, v_min2, v_max2, "yx", out_fldr)
# compareStrain(strains_elastix[1:-1,1:-1,1,1], img_frho3, v_min4, v_max4, "yy", out_fldr)



# compareStrain(strains_elastix, strains_analyt, vmin, vmax, out_fldr)



