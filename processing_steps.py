""" Copyright (C) 2022 E-Scopics. All rights reserved.

    This code is the exclusive property of E-Scopics.
    For research purposes only. Do not redistribute.

    For any question, please contact the author:
    Baptiste Hériard-Dubreuil,
    <baptiste.heriard-dubreuil@e-scopics.com>
"""

import numpy as np
import scipy.interpolate as interp
import scipy.signal as sig

def pw_reception(s_transducers, rx_angle, c_th, pitch, 
                 sampling_frequency, demodulation_frequency,
                 initial_times):
    """
    Synthetizes a plane wave reception from transducer signals,
    based on eq (3).
    Inputs:
        - s_transducers: 2D-array of size MxN corresponding to the 
        signals received by transducers.
        - rx_angle: Desired angle of the received plane-wave (in rad).
        - c_th: Theoretical SOS.
        - pitch: pitch of the probe.
        - sampling_frequency: sampling frequency of s_transducers.
        - demodulation_frequency: demodulation frequency.
        - initial_times: Time instant of the first sample of each 
        transducer signals (w.r.t a time origin chosen at the 
        intersection of the tx pw and the probe center).
    Returns:
        - A 1D-array of size N corresponding to the received plane wave
        of angle rx_angle sampled at sampling_frequency.
    """
    nb_samples = s_transducers.shape[1]
    nb_elements = s_transducers.shape[0]
    time_axis_iq = np.arange(nb_samples) / sampling_frequency
    depth = time_axis_iq / 2 * c_th
    
    # Lateral windowing to select the region [-1cm, 1cm]
    x0 = -1e-2 - depth * np.tan(rx_angle)
    x1 = 1e-2 - depth * np.tan(rx_angle)
    
    s_angular = np.zeros((nb_samples), dtype=np.complex64)
    for i_e in range(nb_elements):
        # Equation (3)
        x_i_e = (i_e - ((nb_elements - 1) / 2)) * pitch
        reception_delay = x_i_e * np.sin(rx_angle) / c_th
        
        # Interpolation
        s_interp_real = np.interp(time_axis_iq - reception_delay, 
                                  time_axis_iq + initial_times[i_e],
                                  np.real(s_transducers[i_e, :]))
        s_interp_imag = np.interp(time_axis_iq - reception_delay, 
                                  time_axis_iq + initial_times[i_e],
                                  np.imag(s_transducers[i_e, :]))
        # Lateral windowing
        s_angular += (s_interp_real + 1j*s_interp_imag) * (x_i_e > x0) * (x_i_e < x1)
        
    return s_angular

def rescale_signal(signal, scaling_factor):
    """
    Rescale signal with a given scaling factor.
    Inputs:
        - signal: A 1D signal.
        - scaling_factor: The wanted scaling facor
    Returns:
        - The re-scaled signal.
    """
    sample_axis = np.arange(len(signal))
    
    s_interp_real = np.interp(sample_axis * scaling_factor, 
                              sample_axis, np.real(signal))
    s_interp_imag = np.interp(sample_axis * scaling_factor, 
                              sample_axis, np.imag(signal))
    
    rescaled_signal = s_interp_real + 1j * s_interp_imag
    
    return rescaled_signal

def measure_scaling_factor(s, s_prime, averaging_length, sampling_frequency):
    """
    Measure the residual scaling factor between two signals using eq (13).
    Inputs:
        - s: The first signal.
        - s_prime: The second signal.
        - averaging_length: The averaging length in number of samples.
        - sampling_frequency: The sampling frequency of s and s_prime.
    Returns:
        - The measured scaling factor.
        - A confidence weight
    """
    # Correlation
    corr = s * s_prime.conj() / (np.abs(s * s_prime) + 1e-16)
    
    # Averaging
    filt = sig.windows.hann(averaging_length)
    corr_avg = sig.convolve(corr, filt, mode='valid')
    
    # Phase
    corr_avg_phase = np.unwrap(np.angle(corr_avg))
    
    # Derivation
    scaling_factor = np.diff(corr_avg_phase) * sampling_frequency
    
    return scaling_factor, np.abs(corr_avg)[:-1]

def sos_from_rsf(rsf, central_frequency, delta_th, delta_th_prime, c_th):
    """
    Compute the local sos from the residual scaling factor according
    to eq (15).
    Inputs:
        - rsf: Residual scaling factor.
        - central_frequency: The sign al central frequency
        - delta_th: The theoretical half-difference of tx and rx angles 
        for the first signal.
        - delta_th_prime: The theoretical half-difference of tx and rx angles 
        for the second signal.
        - c_th: The theoretical sos.
    Returns:
        - The local SOS
    """
    
    sos = c_th * np.sqrt(1 - 2 / (2 * np.pi * central_frequency) * 
                         rsf / (delta_th_prime ** 2 - delta_th ** 2))
    return sos