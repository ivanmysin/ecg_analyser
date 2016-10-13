# -*- coding: utf-8 -*-
"""
ecg analizer
"""
import numpy  as np
import scipy.signal as sig
import matplotlib.pyplot as plt


def get_ecg(capture_length, sampling_rate, noise_varience=0.05, aritmy_noise=0.05, bpm = 60, adc_bit_resolution = 1024):
    # The "Daubechies" wavelet is a rough approximation to a real,
    # single, heart beat ("pqrst") signal
    pqrst = sig.wavelets.daub(10)
    
    # Add the gap after the pqrst when the heart is resting. 
    samples_rest = 10
    zero_array = np.zeros(samples_rest, dtype=float)
    pqrst_full = np.concatenate([pqrst,zero_array])
    
    # Simulated Beats per minute rate
    # For a health, athletic, person, 60 is resting, 180 is intensive exercising

    bps = bpm / 60 # calculate to number bits to seconds
    
    # Simumated period of time in seconds that the ecg is captured in

    
    
    # Caculate the number of beats in capture time period 
    # Round the number to simplify things
    num_heart_beats = int(capture_length * bps)
    
    # Concatonate together the number of heart beats needed
    max_ins = int (sampling_rate * aritmy_noise)
    ecg_template = np.copy(pqrst_full)
    med = np.median(pqrst_full)
    for idx in range(num_heart_beats):
        insN = np.random.randint(0, max_ins) 
        
        ecg_template = np.append(ecg_template, np.random.normal(med, noise_varience, insN) )
        ecg_template = np.append(ecg_template, pqrst_full)
    
    # Add random (gaussian distributed) noise 
    noise = np.random.normal(0, noise_varience, len(ecg_template))
    ecg_template_noisy = noise + ecg_template
    
    # Simulate an ADC by sampling the noisy ecg template to produce the values
    # Might be worth checking nyquist here 
    # e.g. sampling rate >= (2 * template sampling rate)
    num_samples = sampling_rate * capture_length
    ecg_sampled = sig.resample(ecg_template_noisy, num_samples)
    
    # Scale the normalised amplitude of the sampled ecg to whatever the ADC 
    # bit resolution is
    # note: check if this is correct: not sure if there should be negative bit values. 
    
    ecg =  adc_bit_resolution * ecg_sampled
    return ecg

if __name__ == "__main__":
    print('Simulating heart ecg')
    length = 400        # sec
    sampling_rate = 500 # sample in sec
    ecg = get_ecg(length, sampling_rate)  
    
    t = np.linspace(0, length, ecg.size)
    
    # Plot the sampled ecg signal
    plt.figure()
    plt.plot(t, ecg)
    plt.xlim(0, 10)
    plt.xlabel('Time, sec')
    plt.ylabel('Amplitude')
