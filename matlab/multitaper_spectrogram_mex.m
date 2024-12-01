function [mt_spectrogram, stimes, sfreqs] = multitaper_spectrogram_mex(varargin)
%MULTITAPER_SPECTROGRAM_MEX  Compute the multitaper spectrogram for time series data
%
%   Usage:
%   Direct input:
%   [spect,stimes,sfreqs] = multitaper_spectrogram_mex(data, Fs, frequency_range, taper_params, window_params, NFFT, detrend_opt, weighting, plot_on, verbose)
%
%   Input:
%   data: <number of samples> x 1 vector - time series data -- required
%   Fs: double - sampling frequency in Hz  -- required
%
%   (Can be positional or name-value pair arguments)
%   frequency_range: 1x2 vector - [<min frequency>, <max frequency>] (default: [0 nyquist])
%   taper_params: 1x2 vector - [<time-halfbandwidth product>, <number of tapers>] (default: [5 9])
%   window_params: 1x2 vector - [window size (seconds), step size (seconds)] (default: [5 1])
%   NFFT: double - NFFT size, adds zero padding for interpolation (closest 2^x) (default: 2^(nextpow2(<window samples>))
%   detrend_opt: string - detrend data window ('linear' (default), 'constant', 'off');
%   weighting: string - weighting of tapers ('unity' (default), 'eigen', 'adapt');
%   plot_on: boolean to plot results (default: true)
%   verbose: boolean to display spectrogram properties (default: true)
%
%   Output:
%   spect: FxT matrix of spectral power
%   stimes: 1XT vector of times for the center of the spectral bins
%   sfreqs: 1XF vector of frequency bins for the spectrogram
%
%   Example:
%   In this example we create some chirp data and run the multitaper spectrogram on it.
%       Fs=200; %Sampling Frequency
%       frequency_range=[0 25]; %Limit frequencies from 0 to 25 Hz
%       taper_params=[3 5]; %Time bandwidth and number of tapers
%       window_params=[4 1]; %Window size is 4s with step size of 1s
%       nfft=[]; %Use default
%       detrend_opt='constant' %detrend each window by subtracting the average
%       weighting='unity' %weight each taper at 1
%       plot_on=true; %plot spectrogram
%       verbose=true; %print extra info
%
%       %Generate sample chirp data
%       t=1/Fs:1/Fs:600; %Create 10 minutes of data
%       f_start=1;f_end=20; % Set chirp range in Hz
%       data=chirp(t,f_start,t(end),f_end,'logarithmic');
%
%       %Compute the multitaper spectrogram
%       [spect,stimes,sfreqs] = multitaper_spectrogram_mex(data,Fs,frequency_range, taper_params, window_params, nfft, detrend_opt, weighting, plot_on, verbose);
%
%       %Example of name-value pair call
%       [spect,stimes,sfreqs] = multitaper_spectrogram_mex(data, Fs,'nfft', nfft, 'plot_on', true);
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
%% ********************************************************************

%% PROCESS DATA AND PARAMETERS
try
    %Process user input
    [data, Fs, frequency_range, time_bandwidth, num_tapers, winsize_samples, winstep_samples, ~, ~, NFFT, detrend_opt, weighting, plot_on, verbose, xyflip] = process_input(varargin{:});

    if verbose
        display_spectrogram_props([time_bandwidth num_tapers], [winsize_samples winstep_samples], frequency_range, detrend_opt, Fs, NFFT);
    end

    %Generate DPSS tapers (STEP 1)
    [DPSS_tapers, DPSS_eigen] = dpss(winsize_samples, time_bandwidth, num_tapers);

    start_time = datetime('now');
    %Compute the multitaper spectrogram
    if verLessThan('matlab', '9.4') % 9.4 is 2018a
        warning(['Matlab version is ', version('-release'), '. Matlab version must be 2018a or later to run multitaper spectrogram mex. Reverting to matlab version']);
        [mt_spectrogram,stimes,sfreqs] = multitaper_spectrogram(varargin{:});
    else
        [mt_spectrogram, stimes, sfreqs] = multitaper_spectrogram_coder_mex(single(data), Fs, frequency_range, DPSS_tapers, DPSS_eigen, winstep_samples, NFFT, detrend_opt, weighting);
    end

    if xyflip; mt_spectrogram = mt_spectrogram'; end

    %% PLOT THE SPECTROGRAM

    %Show timing if verbose
    if verbose
        disp(' ');
        disp(['Estimation time: ' char(datetime('now')-start_time)]);
    end

    %Plot the spectrogram
    if plot_on
        if xyflip
            imagesc(stimes, sfreqs, nanpow2db(mt_spectrogram'));
        else
            imagesc(stimes, sfreqs, nanpow2db(mt_spectrogram));
        end
        axis xy

        xlabel('Time (s)');
        ylabel('Frequency (Hz)');

        climscale;
        c = colorbar_noresize;
        ylabel(c,'Power (dB)');

        axis tight
    end

catch ME
    warning(ME.message)
    warning(['Mex file execution error for system: ' computer ' . Reverting to matlab version']);
    [mt_spectrogram,stimes,sfreqs] = multitaper_spectrogram(varargin{:});
end
end


% ********************************************
%           HELPER FUNCTIONS
% ********************************************
%% PROCESS THE USER INPUT

function [data, Fs, frequency_range, time_bandwidth, num_tapers, winsize_samples, winstep_samples,...
    window_start, num_windows, nfft, detrend_opt, weighting, plot_on, verbose, xyflip] = process_input(data, Fs, varargin)

% Default values
default_frequency_range = [0 Fs/2];
default_taper_params = [5 9];
default_window_params = [5 1];
default_NFFT_val = 0;
default_detrend_opt = 'linear';
default_weighting = 'unity';
default_plot_on = true;
default_verbose = true;
default_xyflip = false;

% Parse inputs
p = inputParser;
p.addRequired('data', @(x) validateattributes(x, {'numeric','1d'}, {'vector'}));
p.addRequired('Fs', @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
p.addOptional('frequency_range', default_frequency_range, @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x) == 2));
p.addOptional('taper_params', default_taper_params, @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x) == 2));
p.addOptional('window_params', default_window_params, @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x) == 2));
p.addOptional('NFFT_val', default_NFFT_val, @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0 && mod(x, 1) == 0));
p.addOptional('detrend_opt', default_detrend_opt, @(x) isempty(x) || ismember(x, {'linear', 'constant', 'off'}));
p.addOptional('weighting', default_weighting, @(x) isempty(x) || ismember(x, {'unity', 'eigen', 'adapt'}));
p.addOptional('plot_on', default_plot_on, @(x) isempty(x) || (islogical(x) && isscalar(x)));
p.addOptional('verbose', default_verbose, @(x) isempty(x) || (islogical(x) && isscalar(x)));
p.addOptional('xyflip', default_xyflip, @(x) isempty(x) || (islogical(x) && isscalar(x)));

parse(p, data, Fs, varargin{:});

% Access parsed values
data = p.Results.data;
Fs = p.Results.Fs;

% Replace empty values with default values
if isempty(p.Results.frequency_range)
    frequency_range = default_frequency_range;
else
    frequency_range = p.Results.frequency_range;
end

if isempty(p.Results.taper_params)
    taper_params = default_taper_params;
else
    taper_params = p.Results.taper_params;
end

if isempty(p.Results.window_params)
    data_window_params = default_window_params;
else
    data_window_params = p.Results.window_params;
end

if isempty(p.Results.NFFT_val)
    NFFT_val = default_NFFT_val;
else
    NFFT_val = p.Results.NFFT_val;
end

if isempty(p.Results.detrend_opt)
    detrend_opt = 'linear';
else
    detrend_opt = p.Results.detrend_opt;
end

if isempty(p.Results.weighting)
    weighting = 'unity';
else
    weighting = p.Results.weighting;
end

if isempty(p.Results.plot_on)
    plot_on = true;
else
    plot_on = p.Results.plot_on;
end

if isempty(p.Results.verbose)
    verbose = true;
else
    verbose = p.Results.verbose;
end

if isempty(p.Results.xyflip)
    xyflip = false;
else
    xyflip = p.Results.xyflip;
end

%Set either linear or constant detrending
switch lower(detrend_opt)
    case {'const','constant'}
        detrend_opt = 1;
    case 'off'
        detrend_opt = 0;
    otherwise
        detrend_opt = 2;
end

%Set taper weighting options
switch lower(weighting)
    case {'adapt','adaptive'}
        weighting = 2;
    case {'eig', 'eigen'}
        weighting = 1;
    otherwise
        weighting = 0;
end

%Fix error in frequency range
%Set max frequency to nyquist if only lower bound specified
if isscalar(frequency_range)
    frequency_range(2) = Fs/2;
elseif frequency_range(2) > Fs/2
    frequency_range(2) = Fs/2;
    warning(['Upper frequency range greater than Nyquist, setting range to [' num2str(frequency_range(1)) ' ' num2str(frequency_range(2)) ']']);
end

%Set the number of tapers if none supplied
time_bandwidth = taper_params(1);
assert(time_bandwidth>0,'TW must be positive');

%Set the number of tapers to 2 x floor(TW)-1 if none supplied
if isscalar(taper_params)
    num_tapers = floor(2*(time_bandwidth))-1;
    warning(['No taper number specified, setting number of tapers to ' num2str(num_tapers)]);
else
    num_tapers = taper_params(2);
end

%Throw warning for tapers
if num_tapers ~= floor(2*time_bandwidth(1) - 1)
    warning(['Number of tapers is optimal at floor(2*TW - 1). Consider using [' num2str(taper_params(1)) ' ' num2str(floor(2*taper_params(1) - 1)) ']']);
end

%Compute the data window and step size in samples
if mod(data_window_params(1)*Fs,1)
    winsize_samples=round(data_window_params(1)*Fs);
    if mod(data_window_params(1)*Fs,1)>1e-10 %Do not output error if precision issue
        warning(['Window size is not clearly divisible by sampling frequency. Adjusting window size to ' num2str(winsize_samples/Fs) ' seconds']);
    end
else
    winsize_samples=data_window_params(1)*Fs;
end

assert(winsize_samples<=length(data),['Window size is ' num2str(winsize_samples) ' samples but data is only ' num2str(length(data)) ' samples'])
assert(winsize_samples>0, 'Window size must be positive');

if mod(data_window_params(2)*Fs,1)
    winstep_samples=round(data_window_params(2)*Fs);
    if mod(data_window_params(2)*Fs,1)>1e-10 %Do not output error if precision issue
        warning(['Window step size is not clearly divisible by sampling frequency. Adjusting window size to ' num2str(winstep_samples/Fs) ' seconds']);
    end
else
    winstep_samples=data_window_params(2)*Fs;
end

assert(winstep_samples>0, 'Window step must be positive');

%Total data length
N=length(data);

%Force data to be a column vector
if isrow(data)
    data = data(:);
end

%Window start indices
window_start = 1:winstep_samples:N-winsize_samples+1;
%Number of windows
num_windows = length(window_start);

%Number of points in the FFT
if NFFT_val == 0
    nfft = 2^(nextpow2(winsize_samples));
else
    nfft = 2^nextpow2(NFFT_val);
end
end

%% DISPLAY SPECTROGRAM PROPERTIES

function display_spectrogram_props(taper_params, data_window_params, frequency_range, detrend_opt, Fs, NFFT)
data_window_params = data_window_params/Fs;
%my_pool = gcp;
switch detrend_opt
    case 1
        det_string = 'Constant';
    case 2
        det_string = 'Linear';
    otherwise
        det_string='Off';
end

% Display spectrogram properties
disp(' ');
disp('Multitaper Spectrogram Properties:');
disp(' ');
disp(['    Spectral Resolution: ' num2str((2*taper_params(1))/data_window_params(1)) 'Hz']);
disp(['    Window Length: ' num2str(data_window_params(1)) 's']);
disp(['    Window Step: ' num2str(data_window_params(2)) 's']);
disp(['    Time Half-Bandwidth Product: ' num2str(taper_params(1))]);
disp(['    Number of Tapers: ' num2str(taper_params(2))]);
disp(['    Frequency Range: ' num2str(frequency_range(1)) 'Hz - ' num2str(frequency_range(2)) 'Hz']);
disp(['    Frequency Bin Size: ' num2str(Fs/NFFT) 'Hz']);
disp(['    Detrending: ' det_string]);
disp(' ');
end

%% POW2DB FOR NAN ENTRIES

function ydB = nanpow2db(y)
%POW2DB   Power to dB conversion, setting all bad values to nan
%   YDB = POW2DB(Y) convert the data Y into its corresponding dB value YDB
%
%   % Example:
%   %   Calculate ratio of 2000W to 2W in decibels
%
%   y1 = pow2db(2000/2)     % Answer in db

%   Copyright 2006-2014 The MathWorks, Inc.
% EDITED BY MJP 2/7/2020

% #codegen
% cond = all(y(:)>=0);
% if ~cond
%     coder.internal.assert(cond,'signal:pow2db:InvalidInput');
% end

%ydB = 10*log10(y);
%ydB = db(y,'power');
% We want to guarantee that the result is an integer
% if y is a negative power of 10.  To do so, we force
% some rounding of precision by adding 300-300.

ydB = (10.*log10(y)+300)-300;
ydB(y(:)<=0) = nan;
end
