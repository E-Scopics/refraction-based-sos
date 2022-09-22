""" Copyright (C) 2022 E-Scopics. All rights reserved.

    This code is the exclusive property of E-Scopics.
    For research purposes only. Do not redistribute.

    For any question, please contact the author:
    Baptiste Hériard-Dubreuil,
    <baptiste.heriard-dubreuil@e-scopics.com>
"""

from processing_steps import *

def two_signals_sos_estimation(h, tx_angle, h_prime, tx_angle_prime, 
                               gamma, sampling_frequency, demodulation_frequency, 
                               excitation_frequency, c_th, initial_times, 
                               initial_times_prime, pitch, filter_sample_nb):
    # Deduce rx angles and half-differences (delta)
    rx_angle = 2*gamma - tx_angle
    rx_angle_prime = 2*gamma - tx_angle_prime
    delta = gamma - tx_angle
    delta_prime = gamma - tx_angle_prime
    
    # Step 1: PW reception
    s = pw_reception(h, rx_angle, c_th, pitch, sampling_frequency,
                     demodulation_frequency, initial_times)
    s_prime = pw_reception(h_prime, rx_angle_prime, c_th, pitch,
                           sampling_frequency, demodulation_frequency,
                           initial_times_prime)
    
    # Step 2: Rescaling
    s_hat = rescale_signal(s, np.cos(delta))
    s_hat_prime = rescale_signal(s_prime, np.cos(delta_prime))

    # Step 3: RSF Extraction
    rsf, w = measure_scaling_factor(s_hat, s_hat_prime, filter_sample_nb, sampling_frequency)
    
    # Step 4: SOS computation
    sos = sos_from_rsf(rsf, excitation_frequency, delta, delta_prime, c_th)
    
    return sos, w

def multiple_signals_sos_estimation(h_list, tx_angle_list, gamma_list, 
                                    sampling_frequency, demodulation_frequency, 
                                    excitation_frequency, c_th, initial_times_list, 
                                    pitch, filter_sample_nb):
    # Loop over all pairs    
    sos_list = []
    w_list = []
    for i in range(len(h_list)):
        for j in range(i+1, len(h_list)):
            for k in range(len(gamma_list)):
                if check_validity(tx_angle_list[i], tx_angle_list[j], gamma_list[k]):
                    sos, w = two_signals_sos_estimation(h_list[i], tx_angle_list[i], 
                                                     h_list[j], tx_angle_list[j], 
                                                     gamma_list[k], sampling_frequency, 
                                                     demodulation_frequency, 
                                                     excitation_frequency, 
                                                     c_th, initial_times_list[i],
                                                     initial_times_list[j], pitch,
                                                     filter_sample_nb)
                    sos_list.append(sos)
                    w[np.isnan(sos)] = np.nan
                    w_list.append(w)
         
    # Weighted average
    final_sos = np.nansum(np.array(sos_list)*np.array(w_list), axis=0) / np.nansum(w_list, axis=0)
    return final_sos

def check_validity(tx_angle, tx_angle_prime, gamma):
    # Check if the pair is valid (not too high rx angle, delta difference high enough)
    low_rx_angle = np.abs((2*gamma - tx_angle)) < np.deg2rad(30)
    low_rx_prime_angle = np.abs((2*gamma - tx_angle_prime)) < np.deg2rad(30)
    high_delta_diff = np.abs((gamma - tx_angle)**2 - (gamma - tx_angle_prime)**2) > 0.015
    return low_rx_angle * low_rx_prime_angle * high_delta_diff