# -*- coding: utf-8 -*-
"""
analize EGC lib
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelmax
import artifitial_gererator as gen
from itertools import permutations, combinations_with_replacement, combinations, product
from scipy.spatial.distance import mahalanobis


def get_permutations_with_repit(string, N=3):
    arr = []
    for comb in combinations_with_replacement(string, N):
        for perm in permutations(comb, N):
            tmp_str = "" 
            for tmp in perm:
                tmp_str += tmp
            arr.append(tmp_str)
    arr = set(arr)
    return arr

def get_health_by_ecg(ecg, sampling_rate, dataset, verbose=True):
    ecg = np.asarray(ecg).astype(float)
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
        
            
    #print (dcstr)
    codogram = dict()
    codogram_arr = np.array([])
 
    
    codons = get_permutations_with_repit("ABCDEF", 3)

    for codon in codons:
        codogram[codon] = dcstr.count(codon) / (dR.size - 3)
        codogram_arr = np.append(codogram_arr, codogram[codon])

   
    Vinv = np.linalg.inv(np.cov(dataset.T))
    mean_data = np.mean(dataset, axis=0)

    health = mahalanobis(codogram_arr, mean_data, Vinv)
    
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
    dataset = np.loadtxt(datasetfile, delimiter=",") / 597
    
    health = get_health_by_ecg(ecg, sampling_rate, dataset)
    print (health)
    