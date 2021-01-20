# This repository contains 3 multitaper spectrogram implementations (Matlab, R, and Python) from the Prerau Lab

## Matlab Implementation
- multitaper_spectrogram.m: baseline parallelized implementation in Matlab (approximately 5s with default parameters and 8.5hrs of data at 200Hz)
- multitaper_spectrogram_mex.m: optimized implementation in C called from Matlab (approximately 1s with default parameters and 8.5hrs of data at 200Hz)

## Python Implementation
- multitaper_spectrogram_python.py: baseline non-parallelized implementation in Python (approximately 30s with default parameters and 8.5hrs of data at 200Hz)
- requirements.txt: contains names and versions of non-standard library Python packages required to run multitaper_spectrogram_python.py

## R Implementation
- multitaper_spectrogram_R.R: baseline non-parallelized implementation in R (approximately 38s with default parameters and 8.5hrs of data at 200Hz)
