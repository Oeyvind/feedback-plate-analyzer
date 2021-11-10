#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import json
from matplotlib import pyplot as plt
#from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from utils import *
import matplotlib
from matplotlib import cm
cmap = cm.get_cmap("hot")
cmap_r = cm.get_cmap("hot_r")

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
segments = data[6:]
min_logfreq = np.log10(min_freq)
max_logfreq = np.log10(max_freq)

print("Num segments:", len(segments))
#segments=segments[0:14]
#segments=[segments[4]]
#segments=segments[0:5]

columns = min(columns,len(segments))
rows = int(np.ceil(len(segments)/columns))
fig, axes = plt.subplots(rows,1,sharey=True,figsize=(12,16))
segm = 0
for segment in segments:
    col = segm%columns
    row = int(segm/columns)
    if rows > 1: ax = axes[row]
    else: ax = axes
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
    for i in range(len(strongpeaks[0])):
        x=segm+1
        y = strongpeaks[0][i]
        dB_min = -80
        c = cmap((strongpeaks[1][i]/dB_min)+0.1)
        #lineblob(ax,x,y,c,linewidth=1,marker='.')
        ax.hlines(y,x-.4,x+0.4,color=c,linewidth=3)#linewidths=lw)

    #ax.set_xlabel('segm {}'.format(segm))
    if col == 0:
        ax.set_xticks(np.linspace(segm+1,segm+columns,columns))
        ax.set_ylabel('Freq')
        ax.set_ylim(min_freq,max_freq)
        ax.set_yscale('log')
        ax.set_yticks([40,100,300,600,1000,2000,3000])
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2%', pad=0.05)
        cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap_r), cax=cax, orientation='vertical')
        num_ticks = 6
        labels = np.linspace(dB_min, 0, num_ticks)
        cbar.ax.set_yticklabels(["{:d}".format(int(i)) for i in labels]) # add the labels
        cbar.set_label('Amp in dB',labelpad=14, rotation=270)
    segm += 1

plt.tight_layout()
plt.show()


