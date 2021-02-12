# This repository contains 3 multitaper spectrogram implementations (Matlab, R, and Python) from the Prerau Lab

## Matlab Implementation
- **multitaper_spectrogram.m**: baseline parallelized implementation in Matlab (approximately 5s with default parameters and 8.5hrs of data at 200Hz)
- **multitaper_spectrogram_mex.m**: optimized implementation in C called from Matlab (approximately 1s with default parameters and 8.5hrs of data at 200Hz). Data precision is reduced from double to single.

## Python Implementation
- **multitaper_spectrogram_python.py**: implementation in Python with option for parallel processing (approximately 30s with default parameters without parallel processing, approximately 15s with default parameters with parallel processing and 8.5hrs of data at 200Hz)
- **requirements.txt**: contains names and versions of non-standard library Python packages required to run multitaper_spectrogram_python.py
- Results tend to agree with Matlab implementation with precision on the order of at most 10^-12 with SD of at most 10^-10

## R Implementation
- **multitaper_spectrogram_R.R**: baseline non-parallelized implementation in R (approximately 38s with default parameters and 8.5hrs of data at 200Hz)
- Results tend to agree with Python implementation with precision on the order of at most 10^-7 with SD of at most 10^-5
