#!/usr/bin/env python
# -*- coding: utf-8 -*-

# code adapted from https://pythontic.com/visualization/signals/spectrogram
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.io import wavfile
import numpy as np

# Read the wav file (mono)
samplepos = '3_6'
sf = './rec_plate/210816_1446_9x9_explore_pos'+samplepos+'.wav'
sr, signalData = wavfile.read(sf) 
print('sr = ', sr)

fig, ax = plt.subplots(figsize=(10,5))

fftsize = 8192
s,f,t,i=ax.specgram(signalData,Fs=sr, NFFT=fftsize, noverlap=128,norm=colors.PowerNorm(gamma=2), cmap='inferno')
note = ['B','D#','E','F','F#','B','Bb','F','E']
fq = [61,155,164,174,369,246,233,174,164]
xpos = [(2,11),(11,16),(16,20.5),(20.5,26),(26,30),(30,32),(34,38.5),(37,40),(40.5,43)]
for i in range(len(note)):
    ax.plot(xpos[i], (fq[i],fq[i]),color='black')
    if i in [4,6,7]: v='bottom'
    else:v='top'
    ax.text(xpos[i][1]+0.1,fq[i],note[i],color='white', fontweight='bold', verticalalignment=v)
    

ax.axis(ymin=50, ymax=800)
ax.set_yscale('log')
ticks = [60,120,180,240,320,480,720]
ticklabels = [str(ti) for ti in ticks]
ax.set_yticks(ticks)
ax.set_yticklabels(ticklabels)
ax.set_xlabel('Time')
ax.set_ylabel('Frequency')
ax.xaxis.set_minor_locator(MultipleLocator(1))
image_filename = sf[:-3]+'png'
plt.savefig(image_filename)
plt.show()
