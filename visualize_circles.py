#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import json
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from utils import *
import scipy.io.wavfile as wavfile
import sounddevice as sd
import matplotlib
from matplotlib import cm
import colorcet as cc
cmapname="cet_CET_R1_r"
cmap = cm.get_cmap(cmapname)
    
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-sf", "--soundfile", help="sound file name without suffix")
parser.add_argument("-np", "--numpeaks", type=int, help="number of peaks to plot")
parser.add_argument("-v", "--volume", type=int, help="playback volume")
parser.add_argument("-ds", "--distance_semi", type=float, help="minimum peak distance in semitones")
parser.add_argument("-dh", "--distance_Hz", type=float, help="minimum peak distance in Hz")
parser.add_argument("-pn", "--plot_notenames", action="store_true", help="plot pitch class note names of strongest partials")
parser.add_argument("-pc", "--plot_pitchclasses", action="store_true", help="plot pitch class colors of strongest partials")
parser.add_argument("-il", "--plot_interval_lines", action="store_true", help="plot lines showing intervallic structure")
args = parser.parse_args()
if args.soundfile: soundfile = args.soundfile
else: soundfile = '01-in-210601_1556'
if args.numpeaks: numpeaks = args.numpeaks
else: numpeaks = 7
if args.volume: playback_volume = args.volume
else: playback_volume = 0 
plot_notenames = args.plot_notenames
plot_pitchclasses = args.plot_pitchclasses
plot_interval_lines = args.plot_interval_lines

sr, signal = wavfile.read(soundfile+".wav")
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

image_filename = "{}_{}_{}_{}".format(soundfile,cmapname,min_freq,max_freq)
if plot_notenames: image_filename += "_pn"
if plot_pitchclasses: image_filename += "_pc"
if plot_interval_lines: image_filename += "_il"
image_filename += ".png"

print("Num segments:", len(segments))
#segments=segments[0:14]
#segments=[segments[4]]
#segments=segments[0:5]'
segm_strongpeaks = []
segm_fourpeaks = []

columns = min(columns,len(segments))
rows = int(np.ceil(len(segments)/columns))
fig, ax = plt.subplots(figsize=(max(columns,10),max(rows,3)))
segm = 0
notename_labels =[]
pitchclass_squares=[]
interval_lines = []
for segment in segments:
    col = segm%columns+1
    row = int(segm/columns)+1
    segment = peakfilter(segment,distance=distance_semi,distance_type='semitone')
    segment = peakfilter(segment,distance=distance_Hz,distance_type='Hz')
    freqs = np.array(segment[0])
    amps = np.array(segment[1])
    print("segment {} of {}".format(segm+1,len(segments)))
    print("Num peaks: {} for segment {}".format(len(amps),segm+1))
    peaksort = np.argsort(amps)
    strongpeaks = [freqs[peaksort[-numpeaks:]],amps[peaksort[-numpeaks:]]]
    segm_strongpeaks.append(strongpeaks)
    peaksort_freq = np.argsort(strongpeaks[0])
    print("Using the {} strongest peaks".format(numpeaks))
    print("Levels sgm {}".format(segm+1), strongpeaks[0][peaksort_freq],strongpeaks[1][peaksort_freq])
    peakfreq_lognormalized = (np.log10(strongpeaks[0][peaksort_freq])-min_logfreq)/(max_logfreq-min_logfreq)
    print('linewidths', dB_linewidth(strongpeaks[1][peaksort_freq])*0.5)
    cir= circleplot(ax,cmap,col,row,peakfreq_lognormalized,dB_linewidth(strongpeaks[1][peaksort_freq]))
    # strongest 4 peaks
    four_strongpeaks_fq = np.sort(strongpeaks[0][np.argsort(strongpeaks[1])[-min(4,numpeaks):]]) # freq of the four strongest, for display of note names
    segm_fourpeaks.append(four_strongpeaks_fq)
    ## note name labels
    name_offset = 0.36
    labelpositions = [[-name_offset,name_offset],[name_offset,name_offset],[-name_offset,-name_offset],[name_offset,-name_offset]]
    pitchclasses = []
    for pos, fq in enumerate(four_strongpeaks_fq):
        fq_lognormalized = (np.log10(fq)-min_logfreq)/(max_logfreq-min_logfreq)
        labelcolor = cmap(fq_lognormalized)
        xy = np.asarray([col,row])+np.asarray(labelpositions[pos])
        notelabel = ax.annotate(notename(fq),xy=xy, xytext=xy, verticalalignment='center', horizontalalignment='center', color=labelcolor, fontweight='bold')
        notelabel.set_visible(plot_notenames)
        notename_labels.append(notelabel)
        pitchclasses.append(pitchclass(fq))
    ## pitch class color overlays
    overlaycolors = []
    for p in pitchclasses:
        if pitchclasses.count(p) > 1:
            if p not in overlaycolors:
                overlaycolors.append(p)
    print('overlaycolors',overlaycolors)
    if len(overlaycolors) == 1:
        overlay_alpha = 0.5
        overlaycolor = list(cmap(overlaycolors[0]/12))
        overlaycolor[-1] = overlay_alpha
        overlay_rect = patches.Rectangle((col-0.5, row-0.5), 1, 1, linewidth=0, facecolor=overlaycolor)
        overlay_rect.set_visible(plot_pitchclasses)
        pitchclass_squares.append(overlay_rect)
        ax.add_patch(overlay_rect)
    if len(overlaycolors) == 2:
        overlay_alpha1 = 0.5
        overlaycolor1 = list(cmap(overlaycolors[0]/12))
        overlaycolor1[-1] = overlay_alpha1
        overlay_alpha2 = 0.7
        overlaycolor2 = list(cmap(overlaycolors[1]/12))
        overlaycolor2[-1] = overlay_alpha2
        matplotlib.rcParams['hatch.linewidth'] = 4
        overlay_rect = patches.Rectangle((col-0.5, row-0.5), 1, 1, linewidth=0, edgecolor=overlaycolor2, facecolor=overlaycolor1, hatch='//')
        overlay_rect.set_visible(plot_pitchclasses)
        pitchclass_squares.append(overlay_rect)
        ax.add_patch(overlay_rect)
    ## harmonic analysis
    interval_vector = pc_to_intervalvector([pitchclass(fq) for fq in four_strongpeaks_fq], [0,0,0,0,0,0])
    harm_size = 0.3
    harm_pos = intv_xy_pos()
    harm_pos_doubles = 0.03
    for interval in [0,1,2,3,4,5]:
        n_interval = interval_vector[interval]
        i_offsets = intv_offsets(n_interval)
        xoff = [0.5,0.5,1,1,0,0.5]
        yoff = [0.5,0.5,0,0,1,0.5]
        for x,y in zip(np.asarray(i_offsets)*xoff[interval], np.asarray(i_offsets)*yoff[interval]):
            if interval == 5: linestyle=':'
            else: linestyle = '-'
            iline = plt.plot([col+(harm_pos[interval][0][0]*harm_size)+(x*harm_pos_doubles),
                              col+(harm_pos[interval][1][0]*harm_size)+(x*harm_pos_doubles)],
                             [row+(harm_pos[interval][0][1]*harm_size)+(y*harm_pos_doubles),
                              row+(harm_pos[interval][1][1]*harm_size)+(y*harm_pos_doubles)],
                             color='black',linestyle=linestyle)
            iline[0].set_visible(plot_interval_lines)
            interval_lines.append(iline[0])
    segm += 1
    
def onclick(event):
    #print("Click is in pos x {}, y {} ".format(event.xdata, event.ydata))
    sndnum = int(event.xdata-0.5)+int(event.ydata-0.5)*rows
    print('sndnum', sndnum)
    in_corner = False
    #print(np.modf(event.xdata-0.5)[0], np.modf(event.ydata-0.5)[0])
    if (np.modf(event.xdata-0.5)[0]<0.2):
        if (np.modf(event.ydata-0.5)[0]<0.2):
            in_corner = True
            cornerindex = 2
        if (np.modf(event.ydata-0.5)[0]>0.8):
            in_corner = True
            cornerindex = 0
    if (np.modf(event.xdata-0.5)[0]>0.8):
        if (np.modf(event.ydata-0.5)[0]<0.2):
            in_corner = True
            cornerindex = 3
        if (np.modf(event.ydata-0.5)[0]>0.8):
            in_corner = True
            cornerindex = 1
    if in_corner:
        freq = segm_fourpeaks[sndnum][cornerindex]
        print('partial: ', freq, notename(freq))
        soundlen = 2
        t=np.linspace(0,soundlen,sr*soundlen)
        atck = np.linspace(0,1,len(t)//4)
        sust = np.ones(len(t)//2)
        dec = np.linspace(1,0,len(t)//4)
        env = np.concatenate((atck,sust,dec))
        sin =  np.sin(2 * np.pi * freq * t)
        sd.play(sin*0.2*env,sr)
    else:
        if np.modf(event.xdata)[0]<0.5:
            volume = 10**(playback_volume/10)*0.0001
            print('playback volume', playback_volume, 'dB')
            sd.play(signal[audiosegments[sndnum][0]:audiosegments[sndnum][1]]*volume, sr)
        else:
            relative_amp = 1-np.modf(event.ydata-0.5)[0]
            #print('relative amp of partials is scaled by ', relative_amp)
            soundlen = 3
            t=np.linspace(0,soundlen,sr*soundlen)
            s = np.zeros(len(t))
            for i in range(len(segm_strongpeaks[sndnum][0])):
                freq,amp = segm_strongpeaks[sndnum][0][i],segm_strongpeaks[sndnum][1][i]
                atck = np.linspace(0,1,len(t)//4)
                sust = np.ones(len(t)//2)
                dec = np.linspace(1,0,len(t)//4)
                env = np.concatenate((atck,sust,dec))
                sin =  np.sin(2 * np.pi * freq * t)*(10**((amp*relative_amp)/20))*env
                s += sin
                sd.play(s*0.15,sr)
        
cid = fig.canvas.mpl_connect('button_press_event', onclick)

def toggle_notenames(event):
    if event.key not in ['n','p','i']:
        return
    if event.key == 'n':
        for labl in notename_labels:
            labl.set_visible(not labl.get_visible())
        plt.draw()
    if event.key == 'p':
        for pc in pitchclass_squares:
            pc.set_visible(not pc.get_visible())
        plt.draw()
    if event.key == 'i':
        for intv in interval_lines:
            intv.set_visible(not intv.get_visible())
        plt.draw()


plt.connect('key_press_event', toggle_notenames)

divider = make_axes_locatable(ax)
cax = divider.append_axes("left", size="2%", pad=0.4)
norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
mp = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
num_ticks = 10
ticks = np.linspace(0, 1, num_ticks)
labels = np.logspace(np.log10(min_freq), np.log10(max_freq), num_ticks)
cbar = fig.colorbar(mp,cax=cax, orientation='vertical',ticks=ticks)
cax.yaxis.set_ticks_position('left')
cbar.ax.set_yticklabels(["{:d}".format(int(i)) for i in labels]) # add the labels

cax2 = divider.append_axes("right", size="2%", pad=0.2)
dB_list = np.linspace(-70,0,8)
lw = []
for dB in dB_list: lw.append(dB_linewidth(dB))
cax2.hlines(dB_list,-1,1,linewidths=lw)
cax2.yaxis.set_ticks_position('right')
cax2.set_xticks([])

cax3 = divider.append_axes("bottom", size="5%", pad=0.35)
pc_heading = cax3.text(0,1.5,"Pitch class overlay colors")
pitchclass_squares.append(pc_heading)
pc_heading.set_visible(plot_pitchclasses)
pc_alpha1 = 0.5
pc_legend = ['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
for pc in pc_legend:
    pccolor = list(cmap(pc_legend.index(pc)/12))
    pccolor[-1] = pc_alpha1
    pc_rect = patches.Rectangle((pc_legend.index(pc), 0), 0.8, 1, linewidth=0, facecolor=pccolor)
    pitchclass_squares.append(pc_rect)
    pc_rect.set_visible(plot_pitchclasses)
    cax3.add_patch(pc_rect)
    pc_label = cax3.text(0.1+pc_legend.index(pc),0.3,pc)
    pitchclass_squares.append(pc_label)
    pc_label.set_visible(plot_pitchclasses)

intv_heading = cax3.text(12.6,0.7,'Interval legend:\n(Semitones)')
intv_heading.set_visible(plot_interval_lines)
interval_lines.append(intv_heading)
intv_size = 0.4
intv_pos = intv_xy_pos()
intv_x = 16.5
intv_y = 0.5
for interval in [0,1,2,3,4,5]:
    irect = patches.Rectangle((-0.4+intv_x+interval, 0), 0.8, 1, linewidth=0, facecolor=(0.7,0.3,0.2,0.5))
    cax3.add_patch(irect)
    irect.set_visible(plot_interval_lines)
    interval_lines.append(irect)
    if interval == 5: linestyle=':'
    else: linestyle = '-'
    iline = plt.plot([intv_x+interval+(intv_pos[interval][0][0]*intv_size),
                      intv_x+interval+(intv_pos[interval][1][0]*intv_size)],
                     [intv_y+(intv_pos[interval][0][1]*intv_size),
                      intv_y+(intv_pos[interval][1][1]*intv_size)],
                     color='black',linestyle=linestyle)
    iline[0].set_visible(plot_interval_lines)
    interval_lines.append(iline[0])
    intv_label = cax3.text(intv_x+interval,1.3,str(interval+1))
    intv_label.set_visible(plot_interval_lines)
    interval_lines.append(intv_label)
    
cax3.set_xlim([0,22])
cax3.set_ylim([0,2])
cax3.axis('off')

ax.set_aspect(1)
#ax.set_facecolor((0.6,0.6,0.6, 0.2))
ax.set_yticks(np.linspace(1,rows,rows))
ax.set_xticks(np.linspace(1,columns,columns))
plt.tight_layout()
print(image_filename)
plt.savefig(image_filename)
plt.show()


