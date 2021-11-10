#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.colors as colors
from math import log2, pow


def lineblob(p,x,y,c,xwidth=1.0,resolution=100,linewidth=2,marker='o'):
    '''
    Plot a horizontal line that gradially fades to the left and right
    '''
    xs = np.linspace(x-(xwidth*0.5),x+(xwidth*0.5),resolution)
    ys = np.zeros(resolution)
    ys.fill(y)
    a = np.linspace(0,2,resolution)
    a = np.where(a>1, 2-a, a)
    cs = np.empty((resolution,4))
    c = list(c)
    for i in range(len(cs)):
        c[-1] = a[i]
        cs[i] = c
    p.scatter(xs, ys, c=cs, marker=marker, linewidths=linewidth)

def lineblob3D(p,x,y,z,c,xwidth=1.0,resolution=10,linewidth=2,marker='o'):
    '''
    Plot a horizontal line that gradially fades to the left and right
    '''
    xs = np.linspace(x-(xwidth*0.5),x+(xwidth*0.5),resolution)
    ys = np.zeros(resolution)
    ys.fill(y)
    a = np.linspace(0,2,resolution)
    a = np.where(a>1, 2-a, a)
    cs = np.empty((resolution,4))
    c = list(c)
    for i in range(len(cs)):
        c[-1] = a[i]
        cs[i] = c
    p.scatter(xs, ys, z, c=cs, marker=marker, linewidths=linewidth, depthshade=0)
 
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    '''extract part of a colormap'''
    new_cmap = colors.LinearSegmentedColormap.from_list(
       'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
       cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def circleplot(ax,cmap,xorig=0,yorig=0,heights=[0.0,0.5,1.0],widths=[1,2,3]):
    '''
    Plot concentric circles, as horizontal slices of a cone.
    '''
    resolution = 70
    x = np.linspace(-0.5, 0.5, resolution)
    y = np.linspace(-0.5, 0.5, resolution)
    a, b = np.meshgrid(x,y)
    C = np.sqrt(a**2+b**2)
    C = C*(-1/+np.max(C))+1 #invert so center of bump is positive, normalized
    C = C*1.4-0.35 # empirical scaling so heights of 0-1 will look good
    con = ax.contour(a+xorig,b+yorig,C,heights,linewidths=widths,cmap=cmap,norm=colors.Normalize(0,1))
    return con

def circleplot3D(ax,cmap,xorig=0,yorig=0,heights=[0.0,0.5,1.0],widths=[1,2,3]):
    '''
    Plot concentric circles, as horizontal slices of a cone.
    '''
    resolution = 50
    x = np.linspace(-0.5, 0.5, resolution)
    y = np.linspace(-0.5, 0.5, resolution)
    a, b = np.meshgrid(x,y)
    C = np.sqrt(a**2+b**2)
    C = C*(-1/+np.max(C))+1 #invert so center of bump is positive, normalized
    C = C*1.4-0.4 # empirical scaling so heights of 0-1 will look good
    ax.plot_surface(a+xorig,b+yorig,np.clip(C,0,1.1),color=[0,0,0,0.1])
    con = ax.contour(a+xorig,b+yorig,C,heights,linewidths=widths,cmap=cmap,norm=colors.Normalize(0,1))
    #con = ax.plot_surface(a+xorig,b+yorig,np.clip(C,0,1.1),rstride=1, cstride=1, cmap=cmap, linewidth=0, antialiased=False)
    return con

def dB_linewidth(dB):
    ''''
    Convert from dB to linewidth (adjust to taste)
    2**((dB+30)*1/19) gives a linewidth of 0.3 for -60 dB, and width of 3 for -0 dB
    2**((dB+30)*(1/15)) gives a linewidth of 0.25 for -60 dB, and width of 4 for -0 dB
    '''
    return ((2**(dB/12))*2)+1

def dB_ringradius(dB):
    ''''
    Convert from dB to ring radius (adjust to taste)
    -70db should give a radius of 0.003
    0dB should give a radius of 0.03
    2**((dB+30)*1/19) gives a linewidth of 0.3 for -60 dB, and width of 3 for -0 dB
    2**((dB+30)*(1/15)) gives a linewidth of 0.25 for -60 dB, and width of 4 for -0 dB
    '''
    #return 2**((dB-10)*0.04)*0.04
    #return 2**(dB*0.06)*0.05
    return (2**(dB/12)*0.04)+.015

def peakfilter(peaks, distance, distance_type):
    '''
    Filter peaks closer than the specified distance in Hz.
    If two peaks found that are closer, pick the strongest. Check again until no conflicting peaks found.
    
    Parameters
    peaks : list 
        List of two lists, the first containing frequencies, the second containing amplitudes in dB
    distance : int
        The minimum distance allowed between peaks
    distance_type : str
        The scale of the distance value. Can be 'Hz' for linear, or 'semitone' for logarithmic.

    Returns
    outpeaks : list
        List of filtered peaks, same format as input
    '''
    sortindices = np.argsort(peaks[0])
    pops = []
    removed = False
    anyremoved = False
    if distance_type == 'Hz':
        distance_ = distance
    for i in sortindices:
        if i < len(sortindices)-1:
            if not removed:
                if distance_type == 'semitone':
                    distance_ = (peaks[0][i]*np.power(2,distance/12))-peaks[0][i]
                if abs(peaks[0][i]-peaks[0][sortindices[i+1]]) < distance_:
                    removed = True
                    anyremoved = True
                    if peaks[1][i]<peaks[1][sortindices[i+1]]:
                        pops.append(i)
                    else:
                        pops.append(sortindices[i+1])
            else:
                removed = False
    pops.sort()
    pops.reverse()
    for p in pops:
        peaks[0].pop(p)
        peaks[1].pop(p)
    if anyremoved:
        peaks = peakfilter(peaks,distance,distance_type)
    return peaks

def notename(freq):
    '''
    Get the pitch class note name (C,C#,D...) for the given frequency
    
    Parameters
    freq : float
        The frequency to be converted to pitch class note name

    Returns : str
        The note name
    '''
    a4 = 440
    name = ["A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"]
    h = round(12*log2(freq/a4))
    #octave = h // 12
    n = h % 12
    return name[n]

def pitchclass(freq):
    '''
    Get the pitch class normalized of a frequency
    
    Parameters
    freq : float
        The frequency to be converted to pitch class note name

    Returns : int
        The pitch class
    '''
    a4 = 440
    h = round(12*log2(freq/a4))-3
    n = (h % 12)
    return n

def pc_to_intervalvector(pc_list, ivector):
    '''
    Make interval vector from list of pitch classes (pc duplicates are removed)

    Parameters
    pc_list : list
        The list of pitch classes to convert to interval vector
    ivector : list
        The interval vector (under processing).
    Returns : list
        Interval vector
    '''
    pc_list = list(set(pc_list))
    if len(pc_list)>0:
        pc = pc_list.pop(0)
        for i in range(len(pc_list)):
            interval = abs(pc - pc_list[i]) % 12
            if interval != 0:
                if interval > 6: interval = 6 - (interval-6)
                ivector[interval-1] += 1
        pc_to_intervalvector(pc_list, ivector)
    return ivector
    
def intv_offsets(n):
    '''
    Make a list of offsets. If n=2 : [-1,1], if n=3 : [-1,0,1], if n=4[-2,-1,1,2] etc
    
    Parameters
    n : int
        The number of offsets, also the length of the output list.
    Returns : list
        List of offsets, as described above
    '''
    if n%2 == 0:
        offsets = list(np.asarray(range(-int(n/2),int(n/2)))+0.5)
    else:
        offsets = list(range(-int(n/2),int(n/2)+1))
    return offsets

def intv_xy_pos():
    '''
    Make a list of x,y coordinates for plotting intervals as lines intersecting a circle.

    Semitone is a line going from 2 o'clock to 8 o'clock,
    whole tone from 1 to 7, minor third from 12 to 6, etc... tritone from 9 to 3.
    
    Parameters: none
    Returns : list
        List of x and y position offsets
    '''
    '''
    y0 = 0
    y1 = np.sin(np.pi*(1/6))
    y2 = np.sin(np.pi*(2/6))
    y3 = 1
    x0 = 1
    x1 = np.cos(np.pi*(1/6))
    x2 = np.cos(np.pi*(2/6))
    x3 = 0
    xypos = [[[ x1,y1],[-x1,-y1]],
             [[ x2,y2],[-x2,-y2]],
             [[ x3,y3],[-x3,-y3]],
             [[-x2,y2],[ x2,-y2]],
             [[-x1,y1],[ x1,-y1]],
             [[-x0,y0],[ x0,-y0]]]
    '''
    y0 = 0
    y1 = 0.7
    y2 = 1
    x0 = 0
    x1 = 0.7
    x2 = 0.9
    xypos = [[[-x1,-y1],[ x0, y0]],
             [[ x0, y0],[ x1, y1]],
             [[ x0,-y2],[ x0, y0]],
             [[ x0, y0],[ x0, y2]],
             [[-x2, y0],[ x2, y0]],
             [[-x1, y1],[ x1,-y1]]]
    return xypos
