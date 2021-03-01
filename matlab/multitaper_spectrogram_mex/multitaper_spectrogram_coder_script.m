Fs=200; %Sampling Frequency
frequency_range=[0 25]; %Limit frequencies from .5 to 25 Hz
taper_params=[3 5]; %Time bandwidth and number of tapers
window_params=[4 1]; %Window size is 4s with step size of 1s
min_NFFT = 2^10;
detrend_opt = 2;
weighting = 0;

time_bandwidth = taper_params(1);
num_tapers = taper_params(2);

winsize_samples = window_params(1)*Fs;
winstep_samples = window_params(2)*Fs;

%Generate sample chirp data
t=1/Fs:1/Fs:600; %Create 10 minutes of data
f_start=1;f_end=20; % Set chirp range in Hz
data=chirp(t,f_start,t(end),f_end,'logarithmic');


[DPSS_tapers,DPSS_eigen] = dpss(winsize_samples, time_bandwidth, num_tapers);

%Compute the multitaper spectrogram
tic
[spect,stimes,sfreqs] = multitaper_spectrogram_coder(single(data'), Fs, frequency_range, DPSS_tapers, DPSS_eigen, winstep_samples, min_NFFT, detrend_opt, weighting);
toc

figure
imagesc(stimes,sfreqs, pow2db(spect'));
axis xy;
climscale;
colormap(jet);coder
