# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:25:55 2023

@author: laniq
"""
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def secondDeriv(displ_x, displ_y):
    # calculates the second order accuracy derivative using the 
    # definition of derivation
    # returns the second order derivatives in dux/dx, dux/dy, duy/dx, duy/dy using a dictionary
    # where each value in the dictionary is a 2d matrix of derivative values
    
    # in the x displacement field need both the derivative in the x and y direction
    pixelSize = 0.908
    
    # dux/dx gives a 2d matrix
    deriv_x_ydir = (displ_x[1:-1, 2:]-displ_x[1:-1,1:-1])/(2*pixelSize)

    #  dux/dy
    deriv_x_xdir = (displ_x[2:, 1:-1]-displ_x[1:-1,1:-1])/(2*pixelSize)
    
    # in the y displacement field need both the derivative in the x and y direction
    # duy/dx
    deriv_y_ydir = (displ_y[1:-1, 2:]-displ_y[1:-1,1:-1])/(2*pixelSize)
    
    #  duy/dy
    deriv_y_xdir = (displ_y[2:, 1:-1]-displ_y[1:-1,1:-1])/(2*pixelSize)
    
    return {'x_dx': deriv_x_xdir, 'x_dy': deriv_x_ydir, 'y_dx': deriv_y_xdir, 'y_dy': deriv_y_ydir}

def plot_deriv(derivMatrix, outfldr, save):
    fig, axs = plt.subplots(2, 2, sharey=True, figsize=[32,20])
    
    # xx dir
    im1 = axs[0][0].imshow(derivMatrix['x_dx'])
    divider = make_axes_locatable(axs[0][0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
    cbar.set_label('strain', fontsize = 25)
    cbar.ax.tick_params(labelsize=30)
    axs[0][0].axis('off')
    axs[0][0].set_title('Derivative xx direction', fontsize=30)

    # xy dir
    im2 = axs[0][1].imshow(derivMatrix['x_dy'])
    divider = make_axes_locatable(axs[0][1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar1 = fig.colorbar(im2, cax=cax, orientation='vertical')
    cbar1.set_label('strain', fontsize = 25)
    cbar1.ax.tick_params(labelsize=30)
    axs[0][1].axis('off')
    axs[0][1].set_title('Derivative xy direction', fontsize=30)  
    
    # yx dir
    im3 = axs[1][0].imshow(derivMatrix['y_dx'])
    divider = make_axes_locatable(axs[1][0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar2 = fig.colorbar(im3, cax=cax, orientation='vertical')
    cbar2.set_label('strain', fontsize = 25)
    cbar2.ax.tick_params(labelsize=30)
    axs[1][0].axis('off')
    axs[1][0].set_title('Derivative yx direction', fontsize=30)  
    
    # yy dir
    im4 = axs[1][1].imshow(derivMatrix['y_dy'])
    divider = make_axes_locatable(axs[1][1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar3 = fig.colorbar(im4, cax=cax, orientation='vertical')
    cbar3.set_label('strain', fontsize = 25)
    cbar3.ax.tick_params(labelsize=30)
    axs[1][1].axis('off')
    axs[1][1].set_title('Derivative yy direction', fontsize=30)  
    
    plt.suptitle('Estimated Derivatives over Test Image', fontsize = 55)
      
    if save:
        plt.savefig(outfldr + 'derivative_image.png', dpi = 200)
    

def calculate_defGradient(derivMatrix):
    # derivMatrix is a dictionary of second derivatives 
    # where each value in the dictionary is a 2d matrix of derivatives from the original displacement array
    # want to return a matrix of the image where each pixel value is mapped to a matrix of the deformation gradient (f) at that point

    # f for a single pixel is a 2d matrix of four values of the 2nd order derivatives 

    f_matrix = []
    for i in range (len(derivMatrix['x_dx'])):
        row = []
        for j in range (len(derivMatrix['x_dx'][0])):
            sub_matrix = [[(derivMatrix['x_dx'][i][j])+1, derivMatrix['x_dy'][i][j]], [derivMatrix['y_dx'][i][j], (derivMatrix['y_dy'][i][j])+1]]
            row.append(sub_matrix)
        f_matrix.append(np.asarray(row))
        
    return np.asarray(f_matrix)

def calculate_strain(defGradient_matrix):        
    defGradient_matrix = [[0.5*((np.asarray(b).T @ np.asarray(b))-np.identity(2)) for b in a] for a in defGradient_matrix]
        
    return (np.asarray(defGradient_matrix))

def display_strain(strain_matrix, outfldr, save):
    
    fig, axs = plt.subplots(2, 2, sharey=True, figsize=[32,20])
    
    # xx dir
    im1 = axs[0][0].imshow(strain_matrix[:,:,0,0])
    divider = make_axes_locatable(axs[0][0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im1, cax=cax, orientation='vertical');
    cbar.set_label('strain', fontsize = 25)
    cbar.ax.tick_params(labelsize=30)
    axs[0][0].axis('off')
    axs[0][0].set_title('Strain tensor xx direction', fontsize=30)

    # xy dir
    im2 = axs[0][1].imshow(strain_matrix[:,:,0,1])
    divider = make_axes_locatable(axs[0][1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar1 = fig.colorbar(im2, cax=cax, orientation='vertical')
    cbar1.set_label('strain', fontsize = 25)
    cbar1.ax.tick_params(labelsize=30)
    axs[0][1].axis('off')
    axs[0][1].set_title('Strain tensor xy direction', fontsize=30)  
    
    # yx dir
    im3 = axs[1][0].imshow(strain_matrix[:,:,1,0])
    divider = make_axes_locatable(axs[1][0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar2 = fig.colorbar(im3, cax=cax, orientation='vertical')
    cbar2.set_label('strain', fontsize = 25)
    cbar2.ax.tick_params(labelsize=30)
    axs[1][0].axis('off')
    axs[1][0].set_title('Strain tensor yx direction', fontsize=30)  
    
    # yy dir
    im4 = axs[1][1].imshow(strain_matrix[:,:,1,1])
    divider = make_axes_locatable(axs[1][1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar3 = fig.colorbar(im4, cax=cax, orientation='vertical')
    cbar3.set_label('strain', fontsize = 25)
    cbar3.ax.tick_params(labelsize=30)
    # axs[1][1].axis('off')
    axs[1][1].set_title('Strain tensor yy direction', fontsize=30)  
    
    plt.suptitle('Strain Tensor over Test Image', fontsize = 55)
    
    if save:
        plt.savefig(outfldr + 'strain_image.png', dpi = 200)
    
    
    
    
    
    
    
