# -*- coding: utf-8 -*-
"""
analize EGC lib
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelmax
import artifitial_gererator as gen
from itertools import combinations
from scipy.spatial.distance import mahalanobis

def get_health_by_ecg(ecg, sampling_rate, dataset, verbose=True):
    th = 7 * np.median( np.abs(ecg) )
    # out = ecg.ecg(signal=ecg_signal, sampling_rate=1000., show=True)
    
    
    peaks_idx = argrelmax(ecg, order=30)[0]
    peaks_idx = peaks_idx[ ecg[peaks_idx] > th ]
    
    R_amp = ecg[peaks_idx]
    
    T = np.diff(peaks_idx) / sampling_rate
    R_amp = R_amp[0:-1]
    alpha = np.arctan( R_amp / T )
    
    
    
    dR = np.sign (np.diff(R_amp) )    
    dT = np.sign (np.diff(T) )
    dalpha = np.sign (np.diff(alpha) )
    
    dcstr = "" # discret coding string

    for idx in range(dR.size):
        if (dR[idx] > 0 and dT[idx] > 0 and dalpha[idx] > 0):
            dcstr += "A"
            
        if (dR[idx] < 0 and dT[idx] < 0 and dalpha[idx] > 0):
            dcstr += "B"
        
        if (dR[idx] > 0 and dT[idx] < 0 and dalpha[idx] > 0):
            dcstr += "C"
        
        if (dR[idx] < 0 and dT[idx] > 0 and dalpha[idx] < 0):
            dcstr += "D"
        
        if (dR[idx] > 0 and dT[idx] > 0 and dalpha[idx] < 0):
            dcstr += "E"
            
        if (dR[idx] < 0 and dT[idx] < 0 and dalpha[idx] < 0):
            dcstr += "F"
    #dcstr += "ABCDEF"
        
            
    #print (dcstr)
    codogram = dict()
    codogram_arr = np.array([])
    idx = 0
    
    
    for codon in combinations(dcstr, 3):
        key = codon[0] + codon[1] + codon[2]
        codogram[key] = dcstr.count(key) / (dR.size - 3)
        codogram_arr = np.append(codogram_arr, codogram[key])
        print (idx)
        if idx >= 215:
            break
        idx += 1
   
    Vinv = np.linalg.inv(np.cov(dataset.T))
    health = mahalanobis(codogram_arr, np.mean(dataset, axis=0), Vinv)
    
    if (verbose):
        t = np.linspace(0, length, ecg.size)
        
        # Plot the sampled ecg signal
        plt.figure()
        plt.plot(t, ecg)
        plt.plot([t[0], t[-1]], [th, th], "r")
        plt.scatter(t[peaks_idx], ecg[peaks_idx], color="red")
        plt.xlim(0, 10)
        plt.xlabel('Time, sec')
        plt.ylabel('Amplitude')
    
    return (health)


# from biosppy.signals import ecg
if __name__ == "__main__":
    print('Simulating heart ecg')
    length = 600        # sec
    sampling_rate = 500 # sample in sec
    datasetfile = "dataset.csv"
    ecg = gen.get_ecg(length, sampling_rate)
    dataset = np.loadtxt(datasetfile, delimiter=",") / 99
    
    health = get_health_by_ecg(ecg, sampling_rate, dataset)
    print (health)
    