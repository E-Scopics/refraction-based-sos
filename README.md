# Refraction-Based Speed of Sound Estimation in Layered Media

We propose here a code sample that implements the method proposed in:

> B. Hériard-Dubreuil, A. Besson, F. Wintzenrieth, C. Cohen-Bacrie and J.-P. Thiran, "Refraction-Based Speed of Sound Estimation in Layered Media: an Angular Approach" [in preparation].

This code is for research purposes only. No commercial use. Do not redistribute.
If you use the code or the algorithm for research purpose, please cite the above paper.

This code is an explicit version of the code used to obtain the results presented in the above paper. It reproduces the simulation results. For readability reasons, it is not optimized and may take a few minutes to run.

For further data, results reproduction or any question, please contact the author Baptiste Hériard-Dubreuil <baptiste.heriard-dubreuil@e-scopics.com>.


## Code Layout

The proposed code sample is divided in three parts:

- [processing_steps.py](processing_steps.py) that implements the four elementary processing steps described in the paper (see section II.C).

- [sos_estimation.py](sos_estimation.py), with two principal functions, namely **two_signals_sos_estimation** that estimates the SOS from two signals and one mid-angle, as described in section II, and **multiple_signals_sos_estimation** that loops over multiple signals and calls two_signals_sos_estimatio with a lot of different signal pairs and mid-angles (see section II.C).

- [example.ipynb](example.ipynb), a jupyter notebook that uses the above functions and reproduce the paper simulation results.

## Requirements

The following package are required (may work with other versions):

- `numpy==1.18.5`

- `scipy==1.8.1`

- `matplotlib==3.5.2`

To launch **example.ipynb**, please use jupyterlab or jupyter notebook (https://jupyter.org/).


## Copyright

Copyright (C) 2022 E-Scopics. All rights reserved.
This code is the exclusive propriety of E-Scopics.