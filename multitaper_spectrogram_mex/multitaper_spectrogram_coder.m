function [mt_spectrogram,stimes,sfreqs] = multitaper_spectrogram_coder(data, Fs, frequency_range, DPSS_tapers, time_bandwidth, num_tapers, winsize_samples, winstep_samples, min_NFFT, detrend_opt)
%MULTITAPER_SPECTROGRAM  Compute the multitaper spectrogram for time series data
%
%   Usage:
%   Direct input:
%   [spect,stimes,sfreqs] = multitaper_spectrogram_coder(data, Fs, frequency_range, taper_params, window_params, min_NFFT, detrend_opt)
%
%   Input:
%   data: 1 x <number of samples> vector - time series data-- required
%   Fs: double - sampling frequency in Hz  -- required
%   frequency_range: 1x2 vector - [<min frequency>, <max frequency>] (default: [0 nyquist])
%   taper_params: 1x2 vector - [<time-halfbandwidth product>, <number of tapers>] (default: [5 9])
%   window_params: 1x2 vector - [window size (seconds), step size (seconds)] (default: [5 1])
%   detrend_opt: string - detrend data window ('linear' (default), 'constant', 'off');
%   min_NFFT: double - minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0)
%   plot_on: boolean to plot results (default: true)
%   verbose: boolean to display spectrogram properties (default: true)
%
%   Output:
%   spect: TxF matrix of spectral power
%   stimes: 1XT vector of times for the center of the spectral bins
%   sfreqs: 1XF vector of frequency bins for the spectrogram
%
%   Example:
%      Fs=200; %Sampling Frequency
%      frequency_range=[0 25]; %Limit frequencies from .5 to 25 Hz
%      taper_params=[3 5]; %Time bandwidth and number of tapers
%      window_params=[4 1]; %Window size is 4s with step size of 1s
%
%      %Generate sample chirp data
%      t=1/Fs:1/Fs:600; %Create 10 minutes of data
%      f_start=1;f_end=20; % Set chirp range in Hz
%      data=chirp(t,f_start,t(end),f_end,'logarithmic');
%
%      %Compute the multitaper spectrogram
%      [spect,stimes,sfreqs] = multitaper_spectrogram(data,Fs,frequency_range, taper_params, window_params);
%
%   This code is companion to the paper:
%         "Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"
%         Michael J. Prerau, Ritchie E. Brown, Matt T. Bianchi, Jeffrey M. Ellenbogen, Patrick L. Purdon
%         December 7, 2016 : 60-92
%         DOI: 10.1152/physiol.00062.2015
%   which should be cited for academic use of this code.
%
%   A full tutorial on the multitaper spectrogram can be found at:
%   http://www.sleepEEG.org/multitaper
%
%   Copyright 2019 Michael J. Prerau, Ph.D. - http://www.sleepEEG.org
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%   Last modified 1/11/2019
%% ********************************************************************

% PROCESS DATA AND PARAMETERS

%Fix error in frequency range
if frequency_range(2) > Fs/2
    frequency_range(2) = Fs/2;
end

%Generate DPSS tapers (STEP 1)
% DPSS_tapers = coder.const(@dpss, time_bandwidth, num_tapers) * sqrt(Fs);

%Total data length
N=length(data);

%Window start indices
window_start = 1:winstep_samples:N-winsize_samples+1;
%Number of windows
num_windows = length(window_start);

%Number of points in the FFT
nfft = max(max(2^(nextpow2(winsize_samples)),winsize_samples), 2^nextpow2(min_NFFT));

%Create the frequency vector
df = Fs/nfft;
sfreqs = 0:df:Fs; % all possible frequencies

%Set max frequency to nyquist if only lower bound specified
if length(frequency_range) == 1
    frequency_range(2) = Fs/2;
end

%Get just the frequencies for the given frequency range
freq_inds = (sfreqs >= frequency_range(1)) & (sfreqs <= frequency_range(2));
sfreqs = sfreqs(freq_inds);

%Compute the times of the middle of each spectrum
window_middle_times = window_start + round(winsize_samples/2);
stimes = window_middle_times/Fs;

%Preallocate spectrogram and slice data for efficient parallel computing
data_type = class(data);
mt_spectrogram = zeros(sum(freq_inds), num_windows, data_type);
% window_idxs = window_start' + (0:winsize_samples-1);
% data_segments = data(window_idxs);

%% COMPUTE THE MULTITAPER SPECTROGRAM
%
%     STEP 1: Compute DPSS tapers based on desired spectral properties
%     STEP 2: Multiply the data segment by the DPSS Tapers
%     STEP 3: Compute the spectrum for each tapered segment
%     STEP 4: Take the mean of the tapered spectra

%Loop in parallel over all of the windows
parfor n = 1:num_windows
    %Grab the data for the given window
    data_segment = data(window_start(n) + (0:winsize_samples-1));
    
    %Skip empty segments
    if all(data_segment == 0)
        continue;
    end
    
    if any(isnan(data_segment))
        mt_spectrogram(:,n) = nan;
        continue;
    end
    
    %Option to detrend_opt data to remove low frequency DC component
    if detrend_opt == 1
        data_segment = detrend(data_segment, 'constant');
    elseif detrend_opt == 2
        data_segment = detrend(data_segment, 'linear');
    end
    
    %Multiply the data by the tapers (STEP 2)
    tapered_data = repmat(data_segment(:),1,num_tapers) .* DPSS_tapers;
    
    %Compute the FFT (STEP 3)
    fft_data = fft(tapered_data, nfft);
    fft_range = fft_data(freq_inds, :);
    
    %Take the FFT magnitude (STEP 4)
    magnitude = imag(fft_range).^2 + real(fft_range).^2;
    mt_spectrum = sum(magnitude, 2);
    
    %Add the spectrum to the spectrogram
    mt_spectrogram(:,n) = mt_spectrum;
end

%Compute the mean FFT magnitude (STEP 4)
mt_spectrogram = mt_spectrogram' / (Fs*Fs) / num_tapers;

end
