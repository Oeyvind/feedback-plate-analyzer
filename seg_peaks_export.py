#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy.io.wavfile as wavfile
import scipy
import scipy.signal
import numpy as np
from audio_analyze import *
from utils import *
import json
from matplotlib import pyplot as plt
import sounddevice as sd

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-sf", "--soundfile", help="sound file name without suffix")
parser.add_argument("-maxf", "--maxfreq", type=int, help="max freq to plot")
parser.add_argument("-minf", "--minfreq", type=int, help="min freq to plot")
parser.add_argument("-dB", "--dB_thresh", type=int, help="dB threshold for segmentation of audio")
parser.add_argument("-ml", "--min_len", type=float, help="minimum segment length (in seconds)")
parser.add_argument("-rl", "--rms_framelen", type=int, help="The frame length (in milliseconds) over which RMS is measured")
parser.add_argument("-rg", "--relative_gain", action="store_true", help="amps of spectral peaks relative to max amp in file")
parser.add_argument("-ds", "--distance_semi", type=float, help="minimum peak distance in semitones")
parser.add_argument("-dh", "--distance_Hz", type=float, help="minimum peak distance in Hz")
parser.add_argument("-col", "--columns", type=int, help="number of columns and rows (square grid) in display")
parser.add_argument("-d", "--display", help="display/plot data: allowed values are 'wave' or 'spectrum'")
parser.add_argument("-s", "--segments", nargs="+", default=['0','-1'], help="select a subset of segments from the audio input")
parser.add_argument("-ts", "--test_segmentation", action="store_true", help="test segmentation only (no output)")

args = parser.parse_args()
if args.soundfile: soundfile = args.soundfile
else: soundfile = '01-in-210601_1556'
if args.maxfreq: max_freq = args.maxfreq
else: max_freq = 4000
if args.minfreq: min_freq = args.minfreq
else: min_freq = 10
if args.dB_thresh: dB_thresh = args.dB_thresh
else: dB_thresh = -25
if args.min_len: min_len = args.min_len
else: min_len = 1.5
if args.rms_framelen: rms_framelen = args.rms_framelen
else: rms_framelen = 50
if args.relative_gain: use_relative_gain = True
else: use_relative_gain = False
if args.distance_semi: distance_semi = args.distance_semi
else: distance_semi = 1
if args.distance_Hz: distance_Hz = args.distance_Hz
else: distance_Hz = 30
if args.columns: columns = args.columns
else: columns = 9
if args.display: display = args.display
else: display = None
seg_set = args.segments
if args.test_segmentation: test_segmentation = True
else: test_segmentation = False

sr, signal = wavfile.read(soundfile+".wav")
print ("Sampling rate", sr)
l_audio = len(signal.shape)
print ("Num channels", l_audio)
if l_audio == 2:
    print("stereo input, summing to mono")
    signal = signal.sum(axis=1) / 2

N = signal.shape[0]
print ("Number of samples N", N)
secs = N / float(sr)
print ("Sound length in secs", secs)
Ts = 1.0/sr # sampling interval in time
print ("Timestep between samples Ts", Ts)
freq_per_bin = sr/N
print("Frequency resolution", freq_per_bin)
t = np.arange(0, secs, Ts) # time vector as scipy arange field / numpy.ndarray

segments, segments_dB = audio_segmenter(signal,sr,dB_thresh,segment_min_length=min_len, framelength=rms_framelen)
if seg_set != ['0','-1']: 
    segments = segments[int(seg_set[0]):int(seg_set[1])]
print('Num segments:', len(segments))
if test_segmentation:
    print('WARNING: Testing segmentation only (no output)')
    exit()

columns = min(columns,len(segments))
rows = int(np.ceil(len(segments)/columns))
#print('col row', columns,rows)
data = [columns, min_freq, max_freq, dB_thresh, distance_semi,distance_Hz,segments]
if display: fig, axes = plt.subplots(rows,columns,sharey=True,figsize=(max(4,columns*1.5),max(4,rows)))
segm = 0
for segment in segments:
    col = segm%columns
    row = int(segm/columns)
    if display:
        if columns==1:
            ax = axes
        if rows==1:
            ax = axes[col]
        else:
            ax = axes[row,col]
    if use_relative_gain: rel_gain = 10**((np.max(segments_dB) - segments_dB[segm])/20)
    else: rel_gain = 1
    t, fft_freqs, fft_amps_dB, peaks_freqs, peaks_amps_dB = fft_find_peaks(signal[segment[0]:segment[1]], sr, max_freq=max_freq, min_freq=min_freq, p_lookahead=40, p_delta=25,relative_gain=rel_gain)
    #print("segment {} has {} peaks".format(segm+1,len(peaks_freqs)))
    peaks = [peaks_freqs,peaks_amps_dB]
    peaks = peakfilter(peaks,distance=distance_semi,distance_type='semitone')
    peaks = peakfilter(peaks,distance=distance_Hz,distance_type='Hz')
    data.append(peaks)
    if display == 'wave':
        ax.plot(t, signal[segment[0]:segment[1]], "g", alpha=0.7) # plotting the audio segments
    if display == 'spectrum':
        ax.plot(fft_freqs, fft_amps_dB, "g", alpha=0.7) # plotting the spectrum for the audio segments
        ax.scatter(peaks_freqs, peaks_amps_dB, color="r", marker=".") # plotting the spectral peaks for each segment
    if display:
        ax.set_xlabel('{}'.format(segm+1+int(seg_set[0])))
        ax.set_xticks([])
        ax.set_yticks([])

    segm += 1

def flatten_list(_2d_list):
    flat_list = []
    for element in _2d_list:
        try:
            for item in element:
                flat_list.append(item)
        except: flat_list.append(element)
    return flat_list

def onclick(event):
    for i, ax in enumerate(flatten_list(axes)):
        if ax == event.inaxes:
            print("Click is in axes {}".format(i+1))
            soundscale = 4
            sd.play(signal[segments[i][0]:segments[i][1]]*soundscale, sr)

if display: cid = fig.canvas.mpl_connect('button_press_event', onclick)

with open(soundfile+'_peaks.txt', 'w') as outfile:
    json.dump(data, outfile)
'''
with open(soundfile+'_peaks.txt') as json_file:
    data = json.load(json_file)
    for item in data:
        freqs = np.array(item[0])
        amps = np.array(item[1])
        print('freqs', freqs)
'''

if display:
    plt.tight_layout()
    plt.show()
