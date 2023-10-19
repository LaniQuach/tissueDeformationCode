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


trackingpoints = np.load('tracking_points.npy')
xcoords = trackingpoints[:,0,0]
xcoords = [ int(x) for x in xcoords ]
ycoords = trackingpoints[:,0,1]
ycoords = [ int(x) for x in ycoords ]

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
print(I.shape)
ijk = np.vstack([I.flatten(), J.flatten()]).T   # pixel coordinates

xyz_quad, ien_quad, _ = chio.read_mesh(sim_fldr + 'mesh/tissue_quad')
xyz_quad = xyz_quad*1000/pix_vid          # Transform to pixel space

# Read data
frho = chio.read_dfile(sim_fldr + 'data/geomrho.INIT')
srho = chio.read_dfile(sim_fldr + 'data/sarc_rho.INIT')

# Translate mesh
tx = 18/pix_vid
ty = 12/pix_vid
xyz_quad[:,0] += tx
xyz_quad[:,1] += ty

img_frho_y = np.zeros(1502)
img_frho_x = np.zeros(1502)

# Generating a matplotlib triangulation object
tri_quad = Triangulation(xyz_quad[:,0], xyz_quad[:,1], triangles=ien_quad[:,0:3])

numpoints = trackingpoints.shape[0]
tracked0 = np.zeros((numpoints,num_frames+1)) #x dir
tracked1 = np.zeros((numpoints,num_frames+1)) #y dir

for xx in range (0, numpoints):
    tracked0[xx,0] = ycoords[xx]
    tracked1[xx,0] = xcoords[xx]

for ts in tqdm(range(1,num_frames)):
    disp = chio.read_dfile(sim_fldr + 'out/' + 'U-%i' % ts + '.D')
    disp = disp*1000/pix_vid

    interp_func = LinearTriInterpolator(tri_quad, disp[:,0])
    img_frho_x = interp_func(xcoords[:], ycoords[:])
    
    interp_func_y = LinearTriInterpolator(tri_quad, disp[:,1])
    img_frho_y = interp_func_y(xcoords[:], ycoords[:])
    
    print(img_frho_x, img_frho_y)

    for xx in range (0, numpoints): 
        tracked0[xx,ts] = ycoords[xx] + img_frho_y[xx]
        tracked1[xx,ts] = xcoords[xx] + img_frho_x[xx]

# newTrackingPoints = np.zeros([xcoords.shape[0], 1, 2])
# newTrackingPoints[:, 0, 0] = new_points_x
# newTrackingPoints[:, 0, 1] = new_points_y
print(tracked0.shape)
np.savetxt('beat0_row.txt', tracked0)
np.savetxt('beat0_col.txt', tracked1)


# Plot stuff
# plt.figure(0, clear=True)
# plt.imshow(warp_img)
# plt.imshow(img_frho)
# plt.plot(xyz[:,0], xyz[:,1], 'k.')
# plt.plot(fib_xyz[:,0], fib_xyz[:,1], 'k.')
# plt.gca().invert_yaxis()