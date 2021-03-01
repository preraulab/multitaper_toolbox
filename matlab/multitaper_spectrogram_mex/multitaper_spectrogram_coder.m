function [mt_spectrogram, stimes, sfreqs] = multitaper_spectrogram_coder(data, Fs, frequency_range, DPSS_tapers, DPSS_eigen, winstep_samples, min_NFFT, detrend_opt, weighting)
%MULTITAPER_SPECTROGRAM  Compute the multitaper spectrogram for time series data
%
%   This is the coder portion for mex compilation. It takes processed
%   multitaper parameters from multitaper_spectrogram_coder_mex.m as inputs
%   and skips internal input processing. 
%   
%   Usage:
%   Direct input:
%   [spect,stimes,sfreqs] = multitaper_spectrogram_coder(data, Fs, frequency_range, DPSS_tapers, DPSS_eigen, winstep_samples, min_NFFT, detrend_opt, weighting)
%
%   Input:
%   data: <number of samples> x 1 vector - time series data -- required
%   Fs: double - sampling frequency in Hz  -- required
%   frequency_range: 1x2 vector - [<min frequency>, <max frequency>] -- required 
%   DPSS_tapers: Nxk matrix - Slepian tapers -- required 
%   DPSS_eigen: 1xk vector - eigenvalues of the Slepian tapers -- required
%   winstep_samples: double - number of samples in step size of windows -- required 
%   min_NFFT: double - minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) -- required
%   detrend_opt: double - how to detrend data window (2='linear', 1='constant', 0='off') -- required
%   weighting: double - how to weight the tapers (0='unity', 1='eigen', 2='adapt') -- required
%
%   Output:
%   mt_spectrogram: FxT matrix of one-sided power spectral density (PSD)
%   stimes: 1XT vector of times for the center of the spectral bins
%   sfreqs: 1XF vector of frequency bins for the spectrogram
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
%   Last modified 2/16/2021
%% ********************************************************************

%Generate DPSS tapers (STEP 1)
% Done outside this _coder function.

%Get taper matrix dimensions 
[winsize_samples, num_tapers] = size(DPSS_tapers);

%Total data length
N = length(data);

%Window start indices
window_start = 1:winstep_samples:N-winsize_samples+1;

%Number of windows
num_windows = length(window_start);

%Number of points in the FFT
nfft = max(max(2^(nextpow2(winsize_samples)),winsize_samples), 2^nextpow2(min_NFFT));

%Create the frequency vector
df = Fs/nfft;
sfreqs = 0:df:Fs; % all possible frequencies

%Get just the frequencies for the given frequency range
freq_inds = (sfreqs >= frequency_range(1)) & (sfreqs <= frequency_range(2));
sfreqs = sfreqs(freq_inds);

%Compute the times of the middle of each spectrum
window_middle_samples = window_start + round(winsize_samples/2);
stimes = (window_middle_samples-1)/Fs; % stimes start from 0

%Preallocate spectrogram and slice data for efficient parallel computing
data_type = class(data);
mt_spectrogram = zeros(sum(freq_inds), num_windows, data_type);

%% COMPUTE THE MULTITAPER SPECTROGRAM
%
%     STEP 1: Compute DPSS tapers based on desired spectral properties
%     STEP 2: Multiply the data segment by the DPSS Tapers
%     STEP 3: Compute the spectrum for each tapered segment
%     STEP 4: Take the mean of the tapered spectra

% pre-compute weights
if weighting == 1
    wt = DPSS_eigen / num_tapers;
elseif weighting == 0
    wt = ones(num_tapers,1) / num_tapers;
else
    wt = 0;
end

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
    tapered_data = repmat(data_segment,1,num_tapers) .* DPSS_tapers;
    
    %Compute the FFT (STEP 3)
    fft_data = fft(tapered_data, nfft);
    
    %Compute the weighted mean spectral power across tapers (STEP 4)
    Spower = imag(fft_data).^2 + real(fft_data).^2;
    if weighting == 2
        % adaptive weights - for colored noise spectrum (Percival & Walden
        % p368-p370)
        x = data_segment;
        Tpower = x'*x/length(x);
        Spower_iter = mean(Spower(:,1:2),2);
        a = (1-DPSS_eigen)*Tpower;
        for ii = 1:3 % run 3 iterations
            % calculate the MSE weights
            b=(Spower_iter*ones(1,num_tapers))./(Spower_iter*DPSS_eigen'+ones(nfft,1)*a');
            % calculate new spectral estimate
            wk=(b.^2).*(ones(nfft,1)*DPSS_eigen');
            Spower_iter=sum(wk'.*Spower')' ./ sum(wk,2);
        end
        mt_spectrum = Spower_iter;
    else
        % eigenvalue or uniform weights
        mt_spectrum = Spower * wt;
    end
    
    %Append the spectrum to the spectrogram
    mt_spectrogram(:,n) = mt_spectrum(freq_inds);
end

%Compute one-sided PSD spectrum 
DC_select = find(sfreqs==0);
Nyquist_select = find(sfreqs==Fs/2);
select = setdiff(1:length(sfreqs), [DC_select, Nyquist_select]);
mt_spectrogram = [mt_spectrogram(DC_select,:); 2*mt_spectrogram(select,:); mt_spectrogram(Nyquist_select,:)] / Fs;

end
