function [mt_spectrogram,stimes,sfreqs] = multitaper_spectrogram(varargin)
%MULTITAPER_SPECTROGRAM  Compute the multitaper spectrogram for time series data
%
%   Usage:
%   Direct input:
%       [spect,stimes,sfreqs] = multitaper_spectrogram(data, Fs, frequency_range, taper_params, window_params, min_NFFT, detrend_opt, weighting, plot_on, verbose)
%
%   Input:
%       data: <number of samples> x 1  vector - time series data-- required
%       Fs: double - sampling frequency in Hz  -- required
%       frequency_range: 1x2 vector - [<min frequency>, <max frequency>] (default: [0 nyquist])
%       taper_params: 1x2 vector - [<time-halfbandwidth product>, <number of tapers>] (default: [5 9])
%       window_params: 1x2 vector - [window size (seconds), step size (seconds)] (default: [5 1])
%       detrend_opt: string - detrend data window ('linear' (default), 'constant', 'off');
%       min_NFFT: double - minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0)
%       weighting: string - weighting of tapers ('unity' (default), 'eigen', 'adapt');
%       plot_on: boolean to plot results (default: true)
%       verbose: boolean to display spectrogram properties (default: true)
%
%   Output:
%       spect: FxT matrix of spectral power
%       stimes: 1xT vector of times for the center of the spectral bins
%        sfreqs: 1xF vector of frequency bins for the spectrogram
%
%   Example:
%   In this example we create some chirp data and run the multitaper spectrogram on it.
%       Fs=200; %Sampling Frequency
%       frequency_range=[0 25]; %Limit frequencies from 0 to 25 Hz
%       taper_params=[3 5]; %Time bandwidth and number of tapers
%       window_params=[4 1]; %Window size is 4s with step size of 1s
%       min_nfft=0; %No minimum nfft
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
%       [spect,stimes,sfreqs] = multitaper_spectrogram(data,Fs,frequency_range, taper_params, window_params, min_nfft, detrend_opt, weighting, plot_on, verbose);
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
%    Copyright 2021 Michael J. Prerau Laboratory. - http://www.sleepEEG.org
%    Authors: Michael J. Prerau, Ph.D., Mingjian He
%
%   Last modified 1/11/2019
%% ********************************************************************

% PROCESS DATA AND PARAMETERS

%Process user input
[data, Fs, frequency_range, time_bandwidth, num_tapers, winsize_samples, winstep_samples, window_start, num_windows, nfft, detrend_opt, ...
    weighting, plot_on, verbose, xyflip] = process_input(varargin{:});

%Set up and display spectrogram parameters
[window_idxs, stimes, sfreqs, freq_inds] = get_windows(Fs, nfft, frequency_range, window_start, winsize_samples);

if verbose
    display_spectrogram_props([time_bandwidth num_tapers], [winsize_samples winstep_samples], frequency_range, detrend_opt, Fs);
end

%Preallocate spectrogram and slice data for efficient parallel computing
data_type = class(data);
mt_spectrogram = zeros(sum(freq_inds), num_windows, data_type);
data_segments = data(window_idxs)';

%Start timing
start_time = tic;

%% COMPUTE THE MULTITAPER SPECTROGRAM
%
%     STEP 1: Compute DPSS tapers based on desired spectral properties
%     STEP 2: Multiply the data segment by the DPSS Tapers
%     STEP 3: Compute the spectrum for each tapered segment
%     STEP 4: Take the mean of the tapered spectra

%Generate DPSS tapers (STEP 1)
[DPSS_tapers, DPSS_eigen] = dpss(winsize_samples, time_bandwidth, num_tapers);

% pre-compute weights
if weighting == 1
    wt = DPSS_eigen / num_tapers;
elseif weighting == 0
    wt = ones(num_tapers,1) / num_tapers;
else
    wt = 0;
end

%temp_reg = zeros(1024,3,num_windows);

%Loop in parallel over all of the windows
parfor n = 1:num_windows
    %Grab the data for the given window
    data_segment = data_segments(:,n);

    %Skip empty segments
    if all(data_segment == 0)
        continue;
    end

    if any(isnan(data_segment))
        mt_spectrogram(:,n) = nan;
        continue;
    end

    %Option to detrend_opt data to remove low frequency DC component
    if detrend_opt
        data_segment = detrend(data_segment, detrend_opt);
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

    %Add the spectrum to the spectrogram
    mt_spectrogram(:,n) = mt_spectrum(freq_inds);
end


%Compute one-sided PSD spectrum
DC_select = find(sfreqs==0);
Nyquist_select = find(sfreqs==Fs/2);
select = setdiff(1:length(sfreqs), [DC_select, Nyquist_select]);
mt_spectrogram = [mt_spectrogram(DC_select,:); 2*mt_spectrogram(select,:); mt_spectrogram(Nyquist_select,:)] / Fs;

%Flip if requested
if xyflip; mt_spectrogram = mt_spectrogram'; end


%% PLOT THE SPECTROGRAM

%Show timing if verbose
if verbose
    disp(' ');
    disp(['Estimation time: ' datestr(toc(start_time)*datenum([0 0 0 0 0 1]), 'HH:MM:SS.FFF')]);
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

    % Scale color limits of colormap
    climscale;

    % Add colorbar without changing size of figure
    ax = gca;
    pos = ax.Position;
    c = colorbar(ax);
    ax.Position =pos;


    ylabel(c,'Power (dB)');
    axis tight
end

end



% ********************************************
%           HELPER FUNCTIONS
% ********************************************
%% PROCESS THE USER INPUT

function [data, Fs, frequency_range, time_bandwidth, num_tapers, winsize_samples, winstep_samples, window_start, num_windows, nfft, ...
    detrend_opt, weighting, plot_on, verbose, xyflip] = process_input(varargin)
if length(varargin)<2
    error('Too few inputs. Need at least data and sampling rate');
end

%Set default values for inputs
default={[],[],[0 varargin{2}/2],[5 9], [5 1], 0, 'linear', 'unity', true, true, false};

%Allow the third input to be ploton
if nargin == 3 && islogical(varargin{3})
    default{6} = varargin{3};
    varargin = varargin(1:2);
end

%Handle defaults
inputs = default;
inputs(setdiff(1:length(varargin), find(cellfun(@isempty,varargin)))) = varargin(~cellfun(@isempty,(varargin)));

%Transfer input vector to parameters
[data, Fs, frequency_range, taper_params, data_window_params, min_NFFT, detrend_opt, weighting, plot_on, verbose, xyflip] = deal(inputs{:});

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
if length(frequency_range) == 1 %Set max frequency to nyquist if only lower bound specified
    frequency_range(2) = Fs/2;
elseif frequency_range(2) > Fs/2 % updated on 05/18/2020 to remove floor on (Fs/2)
    frequency_range(2) = Fs/2;
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

%Compute the data window and step size in samples
if mod(data_window_params(1)*Fs,1)
    winsize_samples=round(data_window_params(1)*Fs);
    warning(['Window size is not clearly divisible by sampling frequency. Adjusting window size to ' num2str(winsize_samples/Fs) ' seconds']);
else
    winsize_samples=data_window_params(1)*Fs;
end

if mod(data_window_params(2)*Fs,1)
    winstep_samples=round(data_window_params(2)*Fs);
    warning(['Window step size is not clearly divisible by sampling frequency. Adjusting window size to ' num2str(winstep_samples/Fs) ' seconds']);
else
    winstep_samples=data_window_params(2)*Fs;
end

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
nfft = max(max(2^(nextpow2(winsize_samples)),winsize_samples), 2^nextpow2(min_NFFT));
end

%% PROCESS THE SPECTROGRAM PARAMETERS

function [window_idxs, stimes, sfreqs, freq_inds] = get_windows(Fs, nfft, frequency_range, window_start, datawin_size)
%Create the frequency vector
df = Fs/nfft;
sfreqs = 0:df:Fs; % all possible frequencies

%Get just the frequencies for the given frequency range
freq_inds = (sfreqs >= frequency_range(1)) & (sfreqs <= frequency_range(2));
sfreqs = sfreqs(freq_inds);

%Compute the times of the middle of each spectrum
window_middle_samples = window_start + round(datawin_size/2);
stimes = (window_middle_samples-1)/Fs; % stimes start from 0

%Data windows
window_idxs = window_start' + (0:datawin_size-1);
end

%% DISPLAY SPECTROGRAM PROPERTIES

function display_spectrogram_props(taper_params, data_window_params, frequency_range, detrend_opt, Fs)
data_window_params = data_window_params/Fs;
%my_pool = gcp;
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
disp(['    Spectral Resolution: ' num2str((2*taper_params(1))/data_window_params(1)) 'Hz']);
disp(['    Window Length: ' num2str(data_window_params(1)) 's']);
disp(['    Window Step: ' num2str(data_window_params(2)) 's']);
disp(['    Time Half-Bandwidth Product: ' num2str(taper_params(1))]);
disp(['    Number of Tapers: ' num2str(taper_params(2))]);
disp(['    Frequency Range: ' num2str(frequency_range(1)) 'Hz - ' num2str(frequency_range(2)) 'Hz']);
disp(['    Detrending: ' det_string]);
disp(' ');
%disp(['Estimating multitaper spectrogram on ' num2str(my_pool.NumWorkers) ' workers...']);
end

%% POWER TO dB CONVERSION HELPER
function nan_dB = nanpow2db(y)
% Convert power to dB and turn bad values to nan

% We want to guarantee that the result is an integer
% if y is a negative power of 10.  To do so, we force
% some rounding of precision by adding 300-300.
nan_dB = (10.*log10(y)+300)-300;
nan_dB(y(:)<=0) = nan;
end

%% COLOR LIMIT SCALING HELPER FUNCTION
%CLIMSCALE Rescale the color limits of an image to remove outliers with percentiles
%
%   Usage:
%       clim = climscale(hObj, ptiles, outliers)
%       clim(outliers)
%       clim(ptiles)
%
%   Input:
%       hObj: handle to axis or image object -- required
%       ptiles: 1x2 double - scaling percentiles (default: [5 98])
%       outliers: logical - remove outliers prior to scaling using isoutlier (default: true)
%
%   Output:
%       clims: 1x2 double - scaled caxis limits
%
%   Example:
%      ax = gca;
%      imagesc(peaks(500);
%      climscale;
%
%   Copyright 2021 Michael J. Prerau, Ph.D. - http://www.sleepEEG.org
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%   Last modified 10/23/2020
%% ********************************************************************
function clim = climscale(hObj, ptiles, outliers)
if nargin == 1
    if isa(hObj,'matlab.graphics.primitive.Image') || isa(hObj,'matlab.graphics.axis.Axes')
        ptiles =[5 98];
        outliers = true;
    elseif issorted(hObj) && isnumeric(hObj)
        ptiles = hObj;
        hObj = gca;
        outliers = true;
    elseif islogical(hObj)
        outliers = hObj;
        hObj = gca;
        ptiles =[5 98];
    else
        error('Single input must be object, ptiles, or logical');
    end
else
    %Set default current axis
    if nargin==0 || isempty(hObj)
        hObj=gca;
    end

    %Set default percentiles
    if nargin<2 || isempty(ptiles)
        ptiles=[5 98];
    end

    %Set default percentils
    if nargin<3 || isempty(outliers)
        outliers = true;
    end
end

assert(ishandle(hObj) || isa(hObj,'matlab.graphics.primitive.Image') || isa(hObj,'matlab.graphics.axis.Axes'),['First input must be axis or image handle. Input was ' class(hObj)])
assert(issorted(ptiles) && isnumeric(ptiles), 'Percentiles must be monotically increasing and numeric');
assert(islogical(outliers), 'Outliers must be logical');


%Get color data
if isa(hObj,'matlab.graphics.primitive.Image')
    hIm = hObj;
    hAx = get(hIm, 'parent');
else
    hAx = hObj;
    hIm = findall(hAx,'type','image');
    assert(length(hIm) == 1,'More than one image found in axis. Use specific image handle');
end

%Get color data
data = hIm.CData(:);

%Make sure it is not a flat image
assert(range(data)>0,'Image data are all equal');

%Handle massive images
N = length(data);
if N > 1e9
    warning('Data too large to efficiently compute percentile. Using random sampling.');
    data = data(randi(N, 1, min(100000, N)));
end

%Find poorly formed data
if ~outliers
    bad_inds = isnan(data) | isinf(data);
else %Remove outliers if selected
    bad_inds = isnan(data) | isinf(data) | isoutlier(data);
end

%Compute color limits
clim = prctile(data(~bad_inds), ptiles);

if clim(1) == clim(2)
    return;
end

%Update axis scale
set(hAx,'clim',clim);

end

