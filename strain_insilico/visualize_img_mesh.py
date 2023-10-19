#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:52:09 2023

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread, imsave
import cheartio as chio
import math
from matplotlib.tri import Triangulation, LinearTriInterpolator
from tqdm import tqdm

sim_fldr = 'sim_data/'
img_fldr = 'warped_image/images/'
pix_vid = 0.908

# Read image
ts = 25
num_frames = 51
warp_img = imread(img_fldr + 'warped_%02i' % ts + '.png')
i = np.arange(warp_img.shape[0])
j = np.arange(warp_img.shape[1])
I,J = np.meshgrid(i,j)
ijk = np.vstack([I.flatten(), J.flatten()]).T   # pixel coordinates

xyz_quad, ien_quad, _ = chio.read_mesh(sim_fldr + 'mesh/tissue_quad')
xyz_quad = xyz_quad*1000/pix_vid          # Transform to pixel space

# Read data
trackingpoints = np.load('tracking_points.npy')
xcoords = trackingpoints[:,0,0]
xcoords = [ int(x) for x in xcoords ]
ycoords = trackingpoints[:,0,1]
ycoords = [ int(x) for x in ycoords ]

# Translate mesh
tx = 18/pix_vid
ty = 12/pix_vid
xyz_quad[:,0] += tx
xyz_quad[:,1] += ty

numpoints = trackingpoints.shape[0]
tracked0 = np.zeros((numpoints,2)) #x dir
tracked1 = np.zeros((numpoints,2)) #y dir

img_frho_x = np.zeros(1502)
img_frho_y = np.zeros(1502)

# analytical displacement
disp = chio.read_dfile(sim_fldr + 'out/' + 'U-%i' % ts + '.D')
disp = disp*1000/pix_vid

# Generating a matplotlib triangulation object
tri_quad = Triangulation(xyz_quad[:,0], xyz_quad[:,1], triangles=ien_quad[:,0:3])

interp_func = LinearTriInterpolator(tri_quad, disp[:,0])
img_frho_x = interp_func(xcoords[:], ycoords[:])

interp_func_y = LinearTriInterpolator(tri_quad, disp[:,1])
img_frho_y = interp_func_y(xcoords[:], ycoords[:])

#######################################################################
new_points_x = []
new_points_y = []
new_frho_x = []
new_frho_y = []

for p in range(numpoints):
    # get rid of any points that are nan
    if not math.isnan(img_frho_y[p]) and not math.isnan(img_frho_x[p]):
        # update tracking points
        new_points_x.append(xcoords[p])
        new_points_y.append(ycoords[p])
        # update displacement values
        new_frho_x.append(img_frho_x[p])
        new_frho_y.append(img_frho_y[p])

xcoords = np.asarray(new_points_x)
ycoords = np.asarray(new_points_y)
img_frho_x = np.asarray(new_frho_x)
img_frho_y = np.asarray(new_frho_y)
numpoints = ycoords.shape[0]
newTrackingPoints = np.zeros([xcoords.shape[0], 1, 2])
newTrackingPoints[:, 0, 0] = new_points_x
newTrackingPoints[:, 0, 1] = new_points_y

np.save('tracking_points.npy', newTrackingPoints)

for xx in range (0, numpoints):
        tracked0[xx,0] = ycoords[xx]
        tracked1[xx,0] = xcoords[xx]
        tracked0[xx,1] = ycoords[xx] + img_frho_y[xx]
        tracked1[xx,1] = xcoords[xx] + img_frho_x[xx]
     
        
plt.figure()
plt.imshow(warp_img)
plt.scatter(new_points_x[:], new_points_y[:])       
plt.show()

np.savetxt('xdir_AnalyDisp.txt', tracked0)
np.savetxt('ydir_AnalyDisp.txt', tracked1)
