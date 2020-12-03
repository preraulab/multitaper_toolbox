function [mt_spectrum,sfreqs] = multitaper_spectrum(varargin)
%MULTITAPER_SPECTROGRAM  Compute the multitaper spectrogram for time series data
%
%   Usage:
%   Direct input:
%   [spect,stimes,sfreqs] = multitaper_spectrogram(data, Fs, frequency_range, taper_params, min_NFFT, detrend_opt, plot_on, verbose)
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

%Process user input
[data, Fs, frequency_range, time_bandwidth, num_tapers, nfft, detrend_opt, plot_on, verbose] = process_input(varargin{:});
N = length(data);

%Set up and display spectrogram parameters
[sfreqs, freq_inds] = process_spectrogram_params(Fs, nfft, frequency_range);

if verbose
    display_spectrogram_props([time_bandwidth num_tapers], N, frequency_range, detrend_opt, Fs);
end



%Intitalize parallel computing
if isempty(gcp)
    parpool;
end

%Start timing
start_time = tic;

%% COMPUTE THE MULTITAPER SPECTROGRAM
%
%     STEP 1: Compute DPSS tapers based on desired spectral properties
%     STEP 2: Multiply the data segment by the DPSS Tapers
%     STEP 3: Compute the spectrum for each tapered segment
%     STEP 4: Take the mean of the tapered spectra

%Generate DPSS tapers (STEP 1)
DPSS_tapers = dpss(N, time_bandwidth, num_tapers) * sqrt(Fs);

%Option to detrend_opt data to remove low frequency DC component
if detrend_opt
    data = detrend(data, detrend_opt);
end

%Multiply the data by the tapers (STEP 2)
tapered_data = data(:) .* DPSS_tapers;

%Compute the FFT (STEP 3)
fft_data = fft(tapered_data, nfft);
fft_range = fft_data(freq_inds, :);

%Take the FFT magnitude (STEP 4)
magnitude = imag(fft_range).^2 + real(fft_range).^2;
mt_spectrum = sum(magnitude, 2);


%Compute the mean FFT magnitude (STEP 4)
mt_spectrum = mt_spectrum / (Fs*Fs) / num_tapers;

%% PLOT THE SPECTROGRAM

%Show timing if verbose
if verbose
    disp(' ');
    disp(['Estimation time: ' datestr(toc(start_time)*datenum([0 0 0 0 0 1]), 'HH:MM:SS.FFF')]);
end

%Plot the spectrogram
if plot_on
    plot(sfreqs, pow2db(mt_spectrum));
    axis xy
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    axis tight
end
end



% ********************************************
%           HELPER FUNCTIONS
% ********************************************
%% PROCESS THE USER INPUT

function [data, Fs, frequency_range, time_bandwidth, num_tapers, nfft, detrend_opt, plot_on, verbose] = process_input(varargin)
if length(varargin)<2
    error('Too few inputs. Need at least data and sampling rate');
end

%Set default values for inputs
default={[],[],[0 floor(varargin{2}/2)],[5 9], 0, true, true, true};

%Allow the third input to be ploton
if nargin == 3 && islogical(varargin{3})
    default{6} = varargin{3};
    varargin = varargin(1:2);
end

%Handle defaults
inputs = default;
inputs(setdiff(1:length(varargin), find(cellfun(@isempty,varargin)))) = varargin(~cellfun(@isempty,(varargin)));

%Transfer input vector to parameters
[data, Fs, frequency_range, taper_params, min_NFFT, detrend_opt, plot_on, verbose] = deal(inputs{:});

%Set either linear or constant detrending
if detrend_opt ~= false
    switch lower(detrend_opt)
        case {'const','constant'}
            detrend_opt = 'constant';
        case {'none', 'off'}
            detrend_opt = false;
        otherwise
            detrend_opt = 'linear';
    end
end

%Fix error in frequency range
if frequency_range(2) > floor(Fs/2)
    frequency_range(2) = floor(Fs/2);
    warning(['Upper frequency range greater than Nyquist, setting range to [' num2str(frequency_range(1)) ' ' num2str(frequency_range(2)) ']']);
end


%Set the number of tapers if none supplied
time_bandwidth = taper_params(1);

%Set the number of tapers to 2 x floor(TW)-1 if none supplied
if length(taper_params) == 1
    num_tapers = floor(2*(time_bandwidth))-1;
    warning(['No taper number specified, setting number of tapers to ' num2str(num_tapers)]);
else
    num_tapers = taper_params(2);
end

%Throw warning for tapers
if num_tapers ~= floor(2*time_bandwidth(1) - 1)
    warning(['Number of tapers is optimal at floor(2*TW - 1). Consider using [' num2str(taper_params(1)) ' ' num2str(floor(2*taper_params(1) - 1)) ']']);
end

N = length(data);
%Number of points in the FFT
nfft = max(max(2^(nextpow2(N)),N), 2^nextpow2(min_NFFT));
end

%% PROCESS THE SPECTROGRAM PARAMETERS

function [sfreqs, freq_inds] = process_spectrogram_params(Fs, nfft, frequency_range)
%Create the frequency vector
df = Fs/nfft;
sfreqs = 0:df:Fs; % all possible frequencies

%Set max frequency to nyquist if only lower bound specified
if length(frequency_range) == 1
    frequency_range(2) = floor(Fs/2);
end

%Get just the frequencies for the given frequency range
freq_inds = (sfreqs >= frequency_range(1)) & (sfreqs <= frequency_range(2));
sfreqs = sfreqs(freq_inds);

end

%% DISPLAY SPECTROGRAM PROPERTIES

function display_spectrogram_props(taper_params, window_samples, frequency_range, detrend_opt, Fs)
window_time = window_samples/Fs;

if detrend_opt
    det_string=lower(detrend_opt);
    det_string(1) = upper(det_string(1));
else
    det_string='Off';
end

% Display spectrogram properties
disp(' ');
disp('Multitaper Spectrogram Properties:');
disp(' ');
disp(['    Spectral Resolution: ' num2str((2*taper_params(1))/window_time) 'Hz']);
disp(['    Window Length: ' num2str(window_time) 's']);
disp(['    Time Half-Bandwidth Product: ' num2str(taper_params(1))]);
disp(['    Number of Tapers: ' num2str(taper_params(2))]);
disp(['    Frequency Range: ' num2str(frequency_range(1)) 'Hz - ' num2str(frequency_range(2)) 'Hz']);
disp(['    Detrending: ' det_string]);
disp(' ');
end
