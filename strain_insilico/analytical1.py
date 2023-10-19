
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread, imsave
import cheartio as chio
from matplotlib.tri import Triangulation, LinearTriInterpolator

sim_fldr = 'sim_data/'
img_fldr = 'warped_image/images/'
pix_vid = 0.908

# Read image
ts = 45
warp_img = imread(img_fldr + 'warped_%02i' % ts + '.png')
i = np.arange(warp_img.shape[0])
j = np.arange(warp_img.shape[1])
I,J = np.meshgrid(i,j)
ijk = np.vstack([I.flatten(), J.flatten()]).T   # pixel coordinates

# Read mesh
xyz_quad, ien_quad, _ = chio.read_mesh(sim_fldr + 'mesh/tissue_quad')
xyz_quad = xyz_quad*1000/pix_vid          # Transform to pixel space

# Read data
frho = chio.read_dfile(sim_fldr + 'data/geomrho.INIT')
srho = chio.read_dfile(sim_fldr + 'data/sarc_rho.INIT')
C = chio.read_dfile(sim_fldr + 'out/RCauchyGreen-%i' % ts + '.D')
disp = chio.read_dfile(sim_fldr + 'out/' + 'U-%i' % ts + '.D')
disp = disp*1000/pix_vid

# Translate mesh
tx = 18/pix_vid
ty = 12/pix_vid
xyz_quad[:,0] += tx
xyz_quad[:,1] += ty

# Generating a matplotlib triangulation object
tri_quad = Triangulation(xyz_quad[:,0], xyz_quad[:,1], triangles=ien_quad[:,0:3])
interp_func = LinearTriInterpolator(tri_quad, disp[:,0])
img_frho = interp_func(ijk[:,1], ijk[:,0]).data.reshape(I.shape).T

# Plot stuff
plt.figure(0, clear=True)
plt.imshow(warp_img)
plt.imshow(img_frho)
# plt.plot(xyz[:,0], xyz[:,1], 'k.')
# plt.plot(fib_xyz[:,0], fib_xyz[:,1], 'k.')
# plt.gca().invert_yaxis()