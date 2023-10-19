# -*- coding: utf-8 -*-
from skimage.io import imread
from skimage.io import imsave
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


def background():
    # background = np.zeros(400, 600)
    logo_dir = 'data/logo.png'
    
    logo = imread(logo_dir)
    # plt.figure()
    
    xs = np.linspace(-10, 10)
    ys = np.linspace(-10, 10)
    X,Y = np.meshgrid(xs, ys)
    zs = [[np.cos(x) * np.sin(y) for x in xs] for y in ys]
    
    
    background = plt.contourf(X,Y,zs)
    # plt.imshow(logo[:,:])
    plt.axis('off')
    
    # Colorbar(fig[1, 2], hm)
    plt.savefig('data/background.png', dpi = 250)
    
def transpose():
    logo_dir = 'data/logo.png'
    background_dir = 'data/background.png'
    
    logo = Image.open(logo_dir)
    
    image = Image.open(background_dir)
    new_image = image.resize((4000, 3000)) 
    new_image.save("data/background.png")
    extent1 = 0, 4000, 0, 3000
    
    plt.figure()
    plt.imshow(np.asarray(new_image))
    plt.imshow(logo, extent = extent1)

# background()
# 

logo = Image.open('data/targetImage2.png')
new_image = logo.resize((700,500))
new_image.show()
new_image.save('data/targetImage2.png')





    
