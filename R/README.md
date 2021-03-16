# R Multitaper Spectrogram Code
---

<br/>

## Table of Contents
* [General Information](#general-information)
* [Usage](#usage)
* [Example](#example)
* [Parallel Processing](#parallel-processing)
* [Citations](#citations)
* [Status](#status)
* [References](#references)
* [Contact](#contact)

<br/>
<br/>

## General Information 
This folder contains the R implementations of the multitaper spectrogram analysis described in the paper ["Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"](https://prerau.bwh.harvard.edu/publications/Physiology_Bethesda_2017_Prerau.pdf)<sup>1</sup>. Multitaper spectral estimation, which was developed in the early 1980s by David Thomson<sup>2</sup> and has been shown to have superior statistical properties compared with single-taper spectral estimates<sup>3,4</sup>. The multitaper method works by averaging together multiple independent spectra estimated from a single segment of data. The innovation of the multitaper method is that, instead of using a single-taper function to compute the spectrum, it uses multiple taper functions called discrete prolate spheroidal sequences (DPSS). Because DPSS tapers are uncorrelated with each other, they can be averaged together as if they were independent trials of the same condition, producing a spectrum with reduced variance compared to periodogram and single-taper estimation. 

Find videos describing the theory of spectral estimation and demonstrating how multitaper spectral estimation works [http://sleepeeg.org/multitaper](http://sleepeeg.org/multitaper) on the Prerau Lab website. 

<br/>

![alt text](https://prerau.bwh.harvard.edu/images/multitaper_diagram.png)

<sup><sub>Prerau MJ, Bianchi MT, Brown RE, Ellenbogen JM, Patrick PL. Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis. Physiology (Bethesda). 2017 Jan;32(1):60-92. Review. PubMed PMID: 27927806. </sup></sub>

<br/>
<br/>

## Usage
multitaper_spectrogram_R.R utlizies R apply functions, contains no for loops, and includes optional multiprocessing, making it a fairly fast implementation of the multitaper spectrogram analysis.

---

multitaper_spectrogram_R usage:
```
results = multitaper_spectrogram_R(data, fs, frequency_range, time_bandwidth, num_tapers, window_params, min_nfft, weighting, detrend_opt, parallel, num_workers,
                                   plot_on, verbose, xyflip)
spect = results[[1]]
stimes = results[[2]]
sfreqs = results[[3]]
```

<br/>

## Example
In this example we create some chirp data and run the multitaper spectrogram on it.
```
install.packages("signal")
library(signal)  # get signal library for chirp function

# Set spectrogram params
fs = 200  # Sampling Frequency
frequency_range = c(0, 25)  # Limit frequencies from 0 to 25 Hz
time_bandwidth = 3  # Set time-half bandwidth
num_tapers = 5  # Set number of tapers (optimal is time_bandwidth*2 - 1)
window_params = c(4, 1)  # Window size is 4s with step size of 1s
min_nfft = 0  # No minimum nfft
weighting = 'unity'  # weight each taper at 1
detrend_opt = 'constant'  # detrend each window by subtracting the average
parallel = TRUE  # use multiprocessing
num_workers = 3  # use 3 cores in multiprocessing
plot_on = TRUE  # plot spectrogram
verbose = TRUE  # print extra info
xyflip = FALSE  # do not transpose spect output matrix

# Generate sample chirp data
t = seq(1/fs, 600, by=1/fs)  # Create 10 min time array from 1/fs to 600 stepping by 1/fs
f_start = 1  # Set chirp freq range min (Hz)
f_end = 20  # Set chirp freq range max (Hz)
data = chirp(t, f_start, tail(t,n=1), f_end, 'logarithmic')

# Compute the multitaper spectrogram
results = multitaper_spectrogram_R(data, fs, frequency_range, time_bandwidth, num_tapers, window_params, min_nfft, weighting, detrend_opt, parallel, num_workers,
                                   plot_on, verbose, xyflip)
spect = results[[1]]
stimes = results[[2]]
sfreqs = results[[3]]
```
Here is the resulting spectrogram

<img src="https://prerau.bwh.harvard.edu/images/spectrogram_R.png" width="350">

<br/>

## Parallel Processing
The multitaper_spectrogram function makes use of R's 'parallel' and 'doParallel' libraries. To utilize multiprocessing, pass the 'parallel' argument as True and set the 'num_workers' argument to the number of cores you would like to use for parallel processing. If you do not provide the 'num_workers' argument, but the 'parallel' argument is True, the function will automatically use all cores available minus 1. Note that if the 'num_workers' argument exceeds the number of available cores, the function will default to using all available cores minus 1. Also, note that if you choose to use all available cores, your machine will not be able to do anything else while the function is running (because all cores will be in use by the function). Lastly, note that this implementation of parallel processing works fastest on Unix machines (Linux, Mac) due to the way processes are duplicated on these systems. Windows machines will still see benefit from using parallel processing, just not as much as Unix machines.

<br/>

## Citations
The code contained in this repository for multitaper spectral analysis is companion to the paper:  
> "Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"  
>   Michael J. Prerau, Ritchie E. Brown, Matt T. Bianchi, Jeffrey M. Ellenbogen, Patrick L. Purdon  
>    December 7, 2016 : 60-92  
>    DOI: 10.1152/physiol.00062.2015  

which should be cited for academic use of this code.  
<br/>

## Status 
This code is complete and functional but may receive occasional updates. 
<br/>

## References
1. Prerau MJ, Bianchi MT, Brown RE, Ellenbogen JM, Patrick PL. Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis. Physiology (Bethesda). 2017 Jan;32(1):60-92. Review. PubMed PMID: 27927806.
2. Thomson DJ. Spectrum estimation and harmonic analysis. Proc IEEE 70: 1055–1096, 1982.
3. Bronez T. On the performance advantage of multitaper spectral analysis. IEEE Trans Signal Proc 40: 2941–2946, 1992.
4. Percival DB, Walden AT. Spectral Analysis for Physical Applications: Multitaper and Conventional Univariate Techniques. Cambridge, UK: Cambridge Univ. Press, 1993.
<br/>

## Contact
For questions or suggestions please contact Thomas Possidente at tpossidente@bwh.harvard.edu
<br/>

