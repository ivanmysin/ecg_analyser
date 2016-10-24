# -*- coding: utf-8 -*-
"""
ecg classifier by smokers and non-smokers 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelmax


def get_estimate(ecg, sampling_rate, verbose=False):
    ecg = np.asarray(ecg, dtype=float)
    
    th = 1.1 * np.median( np.abs(ecg) )
   
    peaks_idx = argrelmax(ecg, order=5)[0]
    peaks_idx = peaks_idx[ ecg[peaks_idx] > th ]
    
    R_amp = ecg[peaks_idx]
    
    T = np.diff(peaks_idx) / sampling_rate
    #R_amp = R_amp[0:-1]
    heart_rate = np.median (60 / T)
    
    z1 = np.exp ( -(78.40 - heart_rate)**2 / (2 * 11.06**2 ) )# for smokers
    z2 = np.exp( -(71.65 - heart_rate)**2 / (2 * 8.31**2)  )# for non-smokers
    
    if (z1 < 0.1 and z2 < 0.1):
         print ("Classification is not defined")
         print (z1, z2)
         return
    
    p = z1 - z2
    if (p > 0.05 ):
        print ("Smoker")
    elif(p < 0.05 and p > -0.05):
        print ("Classification is not defined")
    elif( p < -0.05):
        print ("Non-smoker")
    print (p)
    
    
    if (verbose):
        t = np.linspace(0, ecg.size/sampling_rate, ecg.size)
        plt.plot(t, ecg, color="b")
        plt.scatter(t[peaks_idx], R_amp, color="red")
        plt.plot([0, t[-1]], [th, th], "r")
        plt.show()

if __name__ == "__main__":
    ecg = np.loadtxt("ecg.txt")
    sampling_rate = 1000
    get_estimate(ecg, sampling_rate, verbose=False)
    

