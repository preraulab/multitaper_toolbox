# Prerau Lab Multitaper Spectrogram Code:
### Matlab, Python, and R implementations
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
* [References](#references)
* [Contact](#contact)

<br/>
<br/>


## Generel Information 
This repository contains Matlab, Python, and R implementations of the multitaper spectrogram analysis described in the paper ["Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"](https://prerau.bwh.harvard.edu/publications/Physiology_Bethesda_2017_Prerau.pdf)<sup>1</sup>. Multitaper spectral estimation, which was developed in the early 1980s by David Thomson<sup>2</sup> and has been shown to have superior statistical properties compared with single-taper spectral estimates<sup>3,4</sup>. The multitaper method works by averaging together multiple independent spectra estimated from a single segment of data. The innovation of the multitaper method is that, instead of using a single-taper function to compute the spectrum, it uses multiple taper functions called discrete prolate spheroidal sequences (DPSS). Because DPSS tapers are uncorrelated with each other, they can be averaged together as if they were independent trials of the same condition, producing a spectrum with reduced variance compared to periodogram and single-taper estimation. 

Find videos describing the theory of spectral estimation and demonstrating how multitaper spectral estimation works [here](https://prerau.bwh.harvard.edu/multitaper/) on the Prerau Lab website. 

<br/>

![alt text](https://prerau.bwh.harvard.edu/images/multitaper_diagram.png)

<sup><sub>Prerau MJ, Bianchi MT, Brown RE, Ellenbogen JM, Patrick PL. Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis. Physiology (Bethesda). 2017 Jan;32(1):60-92. Review. PubMed PMID: 27927806. </sup></sub>

<br/>

## Parameters
The spectral parameters used in all implementations of the multitaper spectrogram are described here
* data: time series data
* fs: sampling frequency in Hz 
* frequency_range: c(<min frequency>, <max frequency>) (default: NULL, adjusted to c(0, nyquist) later)
* taper_params: time_bandwidth: time-half bandwidth product (window duration*half bandwidth of main lobe) (default: 5 Hz*s), num_tapers: number of DPSS tapers to use (default: NULL [will be computed as floor(2*time_bandwidth - 1)])
* window_params: c(window size (seconds), step size (seconds)) (default: [5 1])
* min_nfft (numeric): minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0)
* weighting (char): weighting of tapers ('unity' (default), 'eigen', 'adapt')
* detrend_opt (char): detrend data window ('linear' (default), 'constant', 'off') 

<br/>

## Matlab Implementation
- **multitaper_spectrogram.m**: baseline parallelized implementation in Matlab 
- **multitaper_spectrogram_mex.m**: optimized implementation in C called from Matlab. Data precision is reduced from double to single for major speed improvements.
<br/>

## Python Implementation
- **multitaper_spectrogram_python.py**: baseline implementation in Python with option for multiprocessing
- **requirements.txt**: contains names and versions of non-standard library Python packages required to run multitaper_spectrogram_python.py
- Results tend to agree with Matlab implementation with precision on the order of at most 10^-13 with SD of at most 10^-10
<br/>

## R Implementation
- **multitaper_spectrogram_R.R**: baseline implementation in R (no multiprocessing option implemented yet)
- Results tend to agree with Python implementation with precision on the order of at most 10^-11 with SD of at most 10^-9
<br/>

## Numerical Differences Between Implementations

<br/>

## Weighting Options

<br/>

## Status 
All implementations are complete and functional but may receive updates occasionally
<br/>

## References
1. Prerau MJ, Bianchi MT, Brown RE, Ellenbogen JM, Patrick PL. Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis. Physiology (Bethesda). 2017 Jan;32(1):60-92. Review. PubMed PMID: 27927806.
2. Thomson DJ. Spectrum estimation and harmonic analysis. Proc IEEE 70: 1055–1096, 1982.
3. Bronez T. On the performance advantage of multitaper spectral analysis. IEEE Trans Signal Proc 40: 2941–2946, 1992.
4. Percival DB, Walden AT. Spectral Analysis for Physical Applications: Multitaper and Conventional Univariate Techniques. Cambridge, UK: Cambridge Univ. Press, 1993.
<br/>

## Contact

