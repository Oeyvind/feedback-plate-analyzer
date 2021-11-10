#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy
import numpy as np
from peakdetect import peakdetect

def audio_segmenter(signal, sr, dB_thresh=-20,framelength=50,segment_min_length=1.5):
    '''
    Mark segments in an audio signal, triggered by audio level and minimum segment length.

    Parameters
    signal : array 
        Array of mono audio samples to be analyzed
    sr : int
        The sampling rate of the audio segment
    dB_thresh : int
        The amplitude threshold (in dB), the signal must be louder than this to be included in a segment
    framelength: int
        The RMS frame length in milliseconds
    segment_min_length : float
        The minimum length of a segment in seconds. The signal must be louder than dB_thres for at least segment_min_length

    Returns
    segments : array of lists of ints 
        Segment [start, end] times (indices)
    '''
    N = signal.shape[0]
    segment_start = -1
    segments = [] #indices of each segment
    segments_dB = [] # max amp for each segment
    rms_dB_max = -144
    #signal_rms_dB = np.zeros(len(signal))
    framelen = int(sr*(framelength/1000))
    indx = 0
    while indx < N:
        temp = np.array(signal[indx:indx+framelen],dtype='int64')
        if len(temp) > 0:
            rms = np.sqrt(np.mean(np.square(temp)))
            if rms < 2**-(64) : rms = 2**-(64)
            try:
                rms_dB = 20*np.log10(rms/8192)# not sure why the ref value of 8192 works, but it seems to produce approx correct output
                #print(rms_dB)
                if rms_dB > rms_dB_max: rms_dB_max = rms_dB
            except:
                rms_dB = -144
            if (rms_dB > dB_thresh) and (segment_start < 0):
                segment_start = indx
                rms_dB_max = -144
            if (rms_dB < dB_thresh-6):
                if (segment_start > -1):
                    if ((indx-segment_start) > segment_min_length*sr):
                        segments.append([segment_start, indx])
                        segments_dB.append(rms_dB_max)
                        segment_start = -1
                        rms_dB_max = -144
                else:
                    segment_start = -1
            #rms_dB_a = np.empty(len(temp))
            #rms_dB_a.fill(rms_dB)
        #signal_rms_dB[indx:indx+i_50ms] = rms_dB_a
        indx += framelen

    #print("Num segments :", len(segments))
    return segments, segments_dB

def fft_find_peaks(signal, sr, max_freq=5000, min_freq=20, p_lookahead=150, p_delta=35, relative_gain=1):
    '''
    Find strongest peaks in specrum

    Parameters
    signal : array 
        Array of mono audio samples to be analyzed
    sr : int
        The sampling rate of the audio segment
    max_freq : int
        The maximum frequency to be analyzed and displayed
    relative_gain : float
        Scale the gain by a factor, used to normalize each segment relative to the segment with the max amnplitude
    p_lookahead : int
        The distance (in Hz) to look ahead from a peak to determine if it is the actual peak.
    p_delta : int
        Minimum difference between a peak and the following points

    Returns
    t : array
        Time vector
    fft_freqs : array
        FFT frequencies
    fft_amps_dB : array
        FFT amplitudes in dB
    peaks_freqs : array
        Frequencies of strongest bands
    peaks_amps_dB : array
        Amplitudes (in sdB) of strongest bands
    '''
    N = signal.shape[0]
    secs = N / float(sr)
    Ts = 1.0/sr # sampling interval in time
    freq_per_bin = (sr/N)
    t = np.arange(0, secs, Ts) # time vector as scipy arange field / numpy.ndarray
    if len(t) != N: t=t[:-1]
    w = np.blackman(N)
    FFT = abs(scipy.fft.fft(signal*w))
    FFT_side = FFT[range(N//2)] # one side FFT range
    maxamp = np.amax(abs(FFT_side))*relative_gain
    fft_amps_dB = 20*np.log10(abs(FFT_side)/maxamp)
    fft_freqs = scipy.fft.fftfreq(signal.size, t[1]-t[0])
    fft_freqs = fft_freqs[range(N//2)] # one side frequency range
    fft_freqs = np.array(fft_freqs)
    max_bin = int(max_freq/freq_per_bin)
    min_bin = int(min_freq/freq_per_bin)
    #print('min fq ', fft_freqs[min_bin])
    fft_freqs = fft_freqs[min_bin:max_bin]
    fft_amps_dB = fft_amps_dB[min_bin:max_bin]
    peaks = peakdetect(fft_amps_dB, lookahead=int(p_lookahead/freq_per_bin), delta=p_delta)
    peaks_freqs = [fft_freqs[p[0]] for p in peaks[0]]
    peaks_amps_dB = [fft_amps_dB[p[0]] for p in peaks[0]]
    return t, fft_freqs, fft_amps_dB, peaks_freqs, peaks_amps_dB
