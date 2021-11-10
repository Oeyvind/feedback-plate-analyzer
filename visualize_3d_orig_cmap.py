
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import json
import utils_mlab
from utils import *
from mayavi import mlab
fig = mlab.figure(size=(950,700))
cmap_name = "nipy_spectral"#'gist_rainbow'
c=utils_mlab.color_cone(cmap_name)

cmap = c.module_manager.scalar_lut_manager.lut.table.to_array()
c.module_manager.scalar_lut_manager.reverse_lut = True
shift_index = np.linspace(255*0.07,255*0.93,256)
c.module_manager.scalar_lut_manager.lut.table = cmap[shift_index.astype(int)]

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-sf", "--soundfile", help="sound file name without suffix")
parser.add_argument("-np", "--numpeaks", type=int, help="number of peaks to plot")
parser.add_argument("-ds", "--distance_semi", type=float, help="minimum peak distance in semitones")
parser.add_argument("-dh", "--distance_Hz", type=float, help="minimum peak distance in Hz")
args = parser.parse_args()
if args.soundfile: soundfile = args.soundfile
else: soundfile = '01-in-210601_1556'
if args.numpeaks: numpeaks = args.numpeaks
else: numpeaks = 7

with open(soundfile+'_peaks.txt') as json_file:
    data = json.load(json_file)
columns = data[0]
min_freq = data[1]
max_freq = data[2]
dB_thresh = data[3]
if args.distance_semi: distance_semi = args.distance_semi
else: distance_semi = data[4]
if args.distance_Hz: distance_Hz = args.distance_Hz
else: distance_Hz = data[5]
audiosegments = data[6]
segments = data[7:]
min_logfreq = np.log10(min_freq)
max_logfreq = np.log10(max_freq)

print("Num segments:", len(segments))
#segments=segments[0:14]
#segments=[segments[4]]
#segments=segments[0:5]

columns = min(columns,len(segments))
rows = int(np.ceil(len(segments)/columns))
segm = 0
for segment in segments:
    col = segm%columns
    row = int(segm/columns)
    segment = peakfilter(segment,distance=distance_semi,distance_type='semitone')
    segment = peakfilter(segment,distance=distance_Hz,distance_type='Hz')
    freqs = np.array(segment[0])
    amps = np.array(segment[1])
    print("segment {} of {}".format(segm+1,len(segments)))
    print("Num peaks: {} for segment {}".format(len(amps),segm+1))
    peaksort = np.argsort(amps)
    strongpeaks = [freqs[peaksort[-numpeaks:]],amps[peaksort[-numpeaks:]]]
    #print(strongpeaks)
    peaksort_freq = np.argsort(strongpeaks[0])
    print("Using the {} strongest peaks".format(numpeaks))
    print("Levels sgm {}".format(segm+1), strongpeaks[0][peaksort_freq],strongpeaks[1][peaksort_freq])
    peakfreq_lognormalized = (np.log10(strongpeaks[0][peaksort_freq])-min_logfreq)/(max_logfreq-min_logfreq)
    heights = []
    for i in range(len(peakfreq_lognormalized)):
        heights.append([peakfreq_lognormalized[i], dB_ringradius(strongpeaks[1][peaksort_freq])[i]])
    bgcolor = (0.8,0.8,0.8,0.8)
    s=utils_mlab.plot_cone(col,row,heights,cmap,bgcolor)
    segm += 1

cbar = mlab.colorbar(c, orientation='vertical',label_fmt='')
minpos = 0.145
maxpos = 0.915
def getpos(freq,minpos=minpos,maxpos=maxpos,minfreq=10,maxfreq=4000):
    return ((np.log10(freq)-np.log10(minfreq))/(np.log10(maxfreq)-np.log10(minfreq)))*(maxpos-minpos)+minpos
for f in [10,100,1000,4000]:
    mlab.text(0.055,getpos(f),str(f),width=0.013*len(str(f)))
mlab.view(90, -33, 15, np.array([columns/2-0.7,rows/2-1, 0.0]))
mlab.show()


