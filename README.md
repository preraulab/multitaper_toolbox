# Prerau Lab Multitaper Spectrogram Code
## Matlab, Python, and R implementations
---

## Table of Contents
* [General Information](#gerenal-information)
* [Parameters](#parameters)
* [Matlab Implementation](#matlab-implementation)
* [Python Implementation](#python-implementation)
* [R Implementation](#r-implementation)
* [Numerical Differences Between Implementations](#numerical-differences-between-implementations)
* [Weighting Options](#weighting-options)
* [Status](#status)
* [Contact](#contact)

## Generel Information 

## Parameters

## Matlab Implementation
- **multitaper_spectrogram.m**: baseline parallelized implementation in Matlab 
- **multitaper_spectrogram_mex.m**: optimized implementation in C called from Matlab. Data precision is reduced from double to single.

## Python Implementation
- **multitaper_spectrogram_python.py**: baseline implementation in Python with option for multiprocessing
- **requirements.txt**: contains names and versions of non-standard library Python packages required to run multitaper_spectrogram_python.py
- Results tend to agree with Matlab implementation with precision on the order of at most 10^-12 with SD of at most 10^-10

## R Implementation
- **multitaper_spectrogram_R.R**: baseline implementation in R (no multiprocessing option implemented yet)
- Results tend to agree with Python implementation with precision on the order of at most 10^-7 with SD of at most 10^-5

## Numerical Differences Between Implementations

## Weighting Options

## Status 

## Contact
