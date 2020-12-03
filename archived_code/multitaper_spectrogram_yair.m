function [mt_spectrogram,stimes,sfreqs]=multitaper_spectrogram(data,Fs,varargin)

%MULTITAPER_SPECTROGRAM  Parallel computation of a detrended multitaper spectrogram

%

%   Usage:

%       [spect,stimes,sfreqs]=multitaper_spectrogram(data, Fs, frequency_range, taper_params, window_params, min_NFFT, detrend_on, plot_on)

%

%   Input:

%       data: 1 x <number of samples> vector - time series data-- required

%       Fs: double - sampling frequency in Hz  -- required

%       frequency_range: 1x2 vector - [<min frequency>, <max frequency>] (default: [0 nyquist])

%       taper_params: 1x2 vector - [<time-halfbandwidth product>, <number of tapers>] (default: [5 9])

%       window_params: 1x2 vector - [window size (seconds), step size (seconds)] (default: [5 1])

%       min_NFFT: double - minimum allowable NFFT size, adds zero padding for short windows (closest 2^x) (default: 0)

%       detrend_on: boolean - linear detrend applied to each window (default: true)

%       plot_on: boolean to plot results (default: true)

%

%   Output:

%       spect: FxT matrix of spectral power

%       stimes: 1XT vector of times for the center of the spectral bins

%       sfreqs: 1XF vector of frequency bins for the spectrogram

%

%   Example:

%      Fs=200; %Sampling Frequency

%      frequency_range=[0 25]; %Limit frequencies from .5 to 25 Hz

%      taper_params=[3 5]; %Time bandwidth and number of tapers

%      window_params=[4 1]; %Window size is 4s with step size of 1s

%

%      %Generate sample chirp data

%      t=1/Fs:1/Fs:600; %Create 10 minutes of data

%      f_start=1;f_end=20; % Start at 10Hz, go up to 400Hz

%      data=chirp(t,f_start,t(end),f_end,'logarithmic');

%

%      %Compute the multitaper spectrogram

%      [spect,stimes,sfreqs]=multitaper_spectrogram(data,Fs,frequency_range, taper_params, window_params);

%

%   This code is companion to the paper:

%         "Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"

%         Michael J. Prerau, Ritchie E. Brown, Matt T. Bianchi, Jeffrey M. Ellenbogen, Patrick L. Purdon

%         December 7, 2016 : 60-92

%         DOI: 10.1152/physiol.00062.2015

%   which should be cited for academic use of this code.

%

%   Copyright 2017 Michael J. Prerau, Ph.D., The General Hospital Corporation

%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.

%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)

%

%   Last modified 1/22/2018

%% ********************************************************************



%Check for minimum input

if nargin<2

    error('Too few inputs. Need at least data and sampling rate (Fs)');

end



%Set up parallel computing

if isempty(gcp)

    parpool;

end

p=gcp;



%%  PARSE INPUT

[frequency_range, taper_params, data_window_params, min_NFFT, detrend_on, plot_on]=parse_input(Fs, varargin{:});



%%  SET UP MULTITAPER PARAMETERS



%Compute the data window and step size in samples

if mod(data_window_params(1)*Fs,1)

    datawin_size=round(data_window_params(1)*Fs);

    warning(['Window size is not clearly divisible by sampling frequency. Adjusting window size to ' num2str(datawin_size/Fs) ' seconds']);

else

    datawin_size=data_window_params(1)*Fs;

end



if mod(data_window_params(2)*Fs,1)

    window_step_size=round(data_window_params(2)*Fs);

    warning(['Window step size is not clearly divisible by sampling frequency. Adjusting window size to ' num2str(window_step_size/Fs) ' seconds']);

else

    window_step_size=data_window_params(2)*Fs;

end



%Get the time bandwidth product

time_bandwidth = taper_params(1);

%Set the number of tapers to 2xfloor(TW)-1 if none supplied

if length(taper_params)==1

    num_tapers = floor(2*(time_bandwidth))-1;

else

    num_tapers=taper_params(2);

end



%Total data length

N=length(data);



%Window start indices

window_start=1:window_step_size:N-datawin_size+1;

%Number of windows

num_windows=length(window_start);



%Set the number of points in the FFT

nfft=max(max(2^(nextpow2(datawin_size)),datawin_size), 2^nextpow2(min_NFFT));



%%  DISPLAY SPECTROGRAM PROPERTIES

disp(' ');

disp('Spectrogram Properties:');

disp(' ');

disp(['    Spectral Resolution: ' num2str((2*taper_params(1))/data_window_params(1)) 'Hz']);

disp(['    Window Length: ' num2str(data_window_params(1)) 's']);

disp(['    Window Step: ' num2str(data_window_params(2)) 's']);

disp(['    Time Half-Bandwidth Product: ' num2str(taper_params(1))]);

disp(['    Number of Tapers: ' num2str(taper_params(2))]);

disp(['    Frequency Range: ' strrep(num2str(frequency_range),'  ','-') ' Hz']);

if detrend_on, det_string='On'; else, det_string='Off'; end

disp(['    Linear Detrending: ' det_string]);



%% SET UP SECTROGRAM PARAMETERS



%Create the frequency vector

df=Fs/nfft;

sfreqs=0:df:Fs; % all possible frequencies



%Set max frequency to nyquist if only lower bound specified

if length(frequency_range)==1

    frequency_range(2)=floor(Fs/2);

end



%Get just the frequencies for the given frequency range

freq_inds=sfreqs>=frequency_range(1) & sfreqs<=frequency_range(2);

sfreqs=sfreqs(freq_inds);



%Compute the times of the middle of each spectrum

window_middle_times=window_start+round(datawin_size/2);

stimes=window_middle_times/Fs;

stimes=stimes(1:end-1);



%Preallocate the spectrogram

data_type = class(data);  %Yair: double is more accurate but slower

mt_spectrogram = zeros(sum(freq_inds), num_windows-1, data_type);  %Yair: transposed for speed



%% COMPUTE THE MULTITAPER SPECTROGRAM



%Generate DPSS tapers

DPSS_tapers = dpss(datawin_size,time_bandwidth,num_tapers)*sqrt(Fs);



%Get start time

s_time=tic;



%Display computational time

disp(' ');

disp(['Estimating multitaper spectrogram on ' num2str(p.NumWorkers) ' workers...']);



%Loop in parallel over all of the windows

parfor n=1:num_windows-1

    %Grab the data for the given window

    win_start = window_start(n);  %Yair

    indx = win_start : win_start+datawin_size-1;  %Yair



    %Multiply the detrended data by the tapers

    %{

    if detrend_on

        tapered_data=repmat(detrend(data(indx))',1,num_tapers).*DPSS_tapers;

    else

        tapered_data=repmat(data(indx)',1,num_tapers).*DPSS_tapers;

    end

    %}

    %Yair: replaced the above code (in block comment) with the following:

    data_to_process = data(indx)';

    if all(data_to_process==0)

        continue

    end

    if detrend_on

        data_to_process = detrend(data_to_process);

    end

    tapered_data = data_to_process.*DPSS_tapers; % implicit bsxfun is faster than repmat

    %Yair: end faster code block



    %Compute the FFT

    fft_data = fft(tapered_data,nfft); %Yair: /Fs moved below (outside loop) for speed



    %Take the mean FFT magnitude to get the spectrum for the current window

    fft_data_valid = fft_data(freq_inds,:); %Yair

    %magnitude = abs(fft_data_valid).^2; %Yair was: conj(fft_data).*fft_data

    magnitude = imag(fft_data_valid).^2 + real(fft_data_valid).^2; %Yair

    mt_spectrum = sum(magnitude,2);  % Yair was: mean(magnitude,2); /num_tapers moved below (outside loop)



    %Add to the spectrogram

    mt_spectrogram(:,n) = mt_spectrum; %Yair: row-wise access is *much* faster

end

mt_spectrogram = mt_spectrogram / (Fs*Fs) / num_tapers; %Yair: faster when done en-block, outside loop



%Display computational time

disp(['Estimation time: ' datestr(toc(s_time)*datenum([0 0 0 0 0 1]),'HH:MM:SS')]);



%% PLOT THE SPECTROGRAM



%Plot the spectrogram

if plot_on

    disp('Displaying spectrogram...');

    imagesc(stimes,sfreqs,pow2db(mt_spectrogram));  %Yair: mt_spectrogram is already transposed

    axis xy

    xlabel('Time (s)');

    ylabel('Frequency (Hz)');

    %climscale(gca);

    c=colorbar;

    ylabel(c,'Power (dB)');

    %sec2hms;

    axis tight

end



function [frequency_range, taper_params, data_window_params, min_NFFT, detrend_on, plot_on]=parse_input(Fs,varargin)



%Default Parameters

default={[0 Fs/2],[5 9], [5 1], 0, true, true};



%Handle defaults

inputs=default;

validParamIdxs = ~cellfun(@isempty,varargin); %Yair

inputs(validParamIdxs) = varargin(validParamIdxs); %Yair



%Transfer input vector to parameters

frequency_range    = inputs{1};

taper_params       = inputs{2};

data_window_params = inputs{3};

min_NFFT           = inputs{4};

detrend_on         = inputs{5};

plot_on            = inputs{6};

