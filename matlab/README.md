# Matlab Multitaper Spectrogram Code
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
This folder contains the Matlab implementations of the multitaper spectrogram analysis described in the paper ["Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"](https://prerau.bwh.harvard.edu/publications/Physiology_Bethesda_2017_Prerau.pdf)<sup>1</sup>. Multitaper spectral estimation, which was developed in the early 1980s by David Thomson<sup>2</sup> and has been shown to have superior statistical properties compared with single-taper spectral estimates<sup>3,4</sup>. The multitaper method works by averaging together multiple independent spectra estimated from a single segment of data. The innovation of the multitaper method is that, instead of using a single-taper function to compute the spectrum, it uses multiple taper functions called discrete prolate spheroidal sequences (DPSS). Because DPSS tapers are uncorrelated with each other, they can be averaged together as if they were independent trials of the same condition, producing a spectrum with reduced variance compared to periodogram and single-taper estimation. 

Find videos describing the theory of spectral estimation and demonstrating how multitaper spectral estimation works [http://sleepeeg.org/multitaper](http://sleepeeg.org/multitaper) on the Prerau Lab website. 

<br/>

![alt text](https://prerau.bwh.harvard.edu/images/multitaper_diagram.png)

<sup><sub>Prerau MJ, Bianchi MT, Brown RE, Ellenbogen JM, Patrick PL. Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis. Physiology (Bethesda). 2017 Jan;32(1):60-92. Review. PubMed PMID: 27927806. </sup></sub>

<br/>
<br/>

## Usage
The two functions multitaper_spectrogram and multitaper_spectrogram_mex differ only in speed and data precision. The mex function is implemented in C and is therefor much faster, but reduces data precision from double to single. Both function have identicle inputs.  
<br/>
Additionally, the "mex_files" folder contains the compiled C code necessary to run multitaper_spectrogram_mex and must be on the Matlab path in order to be used. Currently, this code is compiled for 64 bit Mac (.mexmaci64), PC (.mexw64), and Linux (.mexa64) use only. 

---

multitaper_spectrogram usage:
```
[spect,stimes,sfreqs] = multitaper_spectrogram(data, Fs, frequency_range, taper_params, window_params, min_NFFT, detrend_opt, weighting, plot_on, verbose)
```

multitaper_spectrogram_mex usage:
```
[spect,stimes,sfreqs] = multitaper_spectrogram_mex(data, Fs, frequency_range, taper_params, window_params, min_NFFT, detrend_opt, weighting, plot_on, verbose)
```
<br/>

## Example
In this example we create some chirp data and run the multitaper spectrogram on it.
```
Fs=200; %Sampling Frequency
frequency_range=[0 25]; %Limit frequencies from 0 to 25 Hz
taper_params=[3 5]; %Time bandwidth and number of tapers
window_params=[4 1]; %Window size is 4s with step size of 1s
min_nfft=0; %No minimum nfft
detrend_opt='constant'; %detrend each window by subtracting the average
weighting='unity'; %weight each taper at 1
plot_on=true; %plot spectrogram
verbose=true; %print extra info

%Generate sample chirp data
t=1/Fs:1/Fs:600; %Create 10 minutes of data
f_start=1;f_end=20; % Set chirp range in Hz
data=chirp(t,f_start,t(end),f_end,'logarithmic');

%Compute the multitaper spectrogram
[spect,stimes,sfreqs] = multitaper_spectrogram(data,Fs,frequency_range, taper_params, window_params, min_nfft, detrend_opt, weighting, plot_on, verbose);
%Or use multitaper_spectrogram_mex
% [spect,stimes,sfreqs] = multitaper_spectrogram_mex(data,Fs,frequency_range, taper_params, window_params, min_nfft, detrend_opt, weighting, plot_on, verbose);
```
Here is the resulting spectrogram

<img src="https://prerau.bwh.harvard.edu/images/chirp_spectrogram.jpg" width="350">

<br/>

## Parallel Processing
The multitaper_spectrogram function makes use of Matlab's parallel processing functionality. Parallel processing is automatically used with the default number of cores, but some may want to initialize parallel processing before calling the function, or set a specific number of cores to be used instead of the default. The code below shows how to do this.
```
p = parcluster; %Get information about current cluster
num_pools = p.NumWorkers; %Set to max workers, or change to a specific number of workers
parpool(num_pools); %Instantiate workers
```
OR
```
gcp; %Start workers with default number for cluster
```
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
This code is complete and functional but may receive updates occasionally
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
