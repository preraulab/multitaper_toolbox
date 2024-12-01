function [mt_spectrum, sfreqs] = multitaper_spectrum_mex(data, Fs, frequency_range, taper_params, varargin)
%MULTITAPER_SPECTRUM_MEX  Compute the multitaper spectrum for time series data
%
%   Usage:
%   Direct input:
%   [spect,stimes,sfreqs] = multitaper_spectrum_mex(data, Fs, frequency_range, taper_params, min_NFFT, detrend_opt, weighting, plot_on, verbose)
%
%   Input:
%   data: <number of samples> x 1 vector - time series data -- required
%   Fs: double - sampling frequency in Hz  -- required
%   frequency_range: 1x2 vector - [<min frequency>, <max frequency>] (default: [0 nyquist])
%   taper_params: 1x2 vector - [<time-halfbandwidth product>, <number of tapers>] (default: [5 9])
%   min_NFFT: double - minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0)
%   detrend_opt: string - detrend data window ('linear' (default), 'constant', 'off');
%   weighting: string - weighting of tapers ('unity' (default), 'eigen', 'adapt');
%   plot_on: boolean to plot results (default: true)
%   verbose: boolean to display spectrogram properties (default: true)
%
%   Output:
%   spect: FxT matrix of spectral power
%   sfreqs: 1XF vector of frequency bins for the spectrogram
%
%   Example:
%       Fs = 200;
%       t = (1/Fs):(1/Fs):100;
%       y = 2*sin(2*pi*t*15) + sin(2*pi*t*10) + .3*sin(2*pi*t*5);
%       [spect, sfreqs] = multitaper_spectrum_mex(y,Fs,[.5 30],[5 9]);
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
%   Copyright 2024 Michael J. Prerau Laboratory. - http://www.sleepEEG.org
%   Authors: Michael J. Prerau, Ph.D., Mingjian He, Ph.D.
%
%% ********************************************************************
%% PROCESS DATA AND PARAMETERS
assert(nargin>= 2, 'Must include data and sampling rate (Fs)');

if nargin<3
    frequency_range = [];
end

if nargin<4
    taper_params = [];
end

%Set the spectrum to be one big window
window_params = [length(data)/Fs length(data)/Fs];

%Set default for plot_on
if length(varargin)<4 || isempty(varargin{4})
    plot_on = true;
else
    plot_on = varargin{4};
end

%Do not do the spectrogram plotting
varargin{4} = false;

try
    [mt_spectrum, ~, sfreqs] = multitaper_spectrogram_mex(data, Fs, frequency_range, taper_params, window_params, varargin{:});
catch ME
    warning(ME.message)
    warning(['Mex file execution error for system: ' computer ' . Reverting to matlab version']);
    [mt_spectrum, ~, sfreqs] = multitaper_spectrogram(data, Fs, frequency_range, taper_params, window_params, varargin{:});
end

if plot_on
    plot(sfreqs,nanpow2db(mt_spectrum));
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    grid on
    title('Multitaper Power Spectral Density')
end