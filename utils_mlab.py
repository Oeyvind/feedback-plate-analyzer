#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from mayavi import mlab

def color_cone(cmap_name):
    resolution = 0.1
    x, y = np.mgrid[-0.5:0.5+resolution:resolution, -0.5:0.5+resolution:resolution]
    z = np.sqrt(x**2+y**2)
    z = z*(-1/+np.max(z))+1 #invert so center of bump is positive, normalized
    z = np.clip(z*1.4-0.4,0,1) # empirical scaling and clipping to cone 
    r = np.linspace(0,np.pi*2,100)
    n = np.zeros(100)
    c = mlab.surf(x,y,z,colormap=cmap_name)
    c.visible = False
    return c
    
def plot_cone(xorig,yorig,items,cmap,bgcolor=[0.5,0.5,0.5,0.7]):
    resolution = 0.05
    x, y = np.mgrid[-0.5:0.5+resolution:resolution, -0.5:0.5+resolution:resolution]
    z = np.sqrt(x**2+y**2)
    z = z*(-1/+np.max(z))+1 #invert so center of bump is positive, normalized
    z = np.clip(z*1.4-0.4,0,1) # empirical scaling and clipping to cone 
    r = np.linspace(0,np.pi*2,100)
    n = np.zeros(100)

    s = mlab.surf(x+xorig,y+yorig,z,color=bgcolor[0:3],opacity=bgcolor[3])
    for item in items:
        height,width = item
        h = (0.5-height*0.5)-0.01+(1-height)*0.04
        mlab.plot3d(np.sin(r)*h+xorig, np.cos(r)*h+yorig, n+height,tube_radius=width, color=tuple(list(cmap[int(height*255)]/255)[:3]))
    return s

