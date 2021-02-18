% Compare multitaper estimation methods 
close all; clear all;

% set up simualtion data 
%duration = 60; %60s
%Fs = 200; %100Hz
%EEG = single(pinknoise(Fs*duration))';
%nw = 2;

[EEG, Fs, ~, ~, ~] = load_Lunesta_data (1, 1, false, 'O1');
nw = 15;

% params
spectrogram_parameters.frequency_max = Fs/2;
spectrogram_parameters.taper_params = [nw, 2*nw-1];
spectrogram_parameters.window_params = [30, 10];
spectrogram_parameters.min_NFFT = 2^10;
spectrogram_parameters.detrend = 'constant';
spectrogram_parameters.ploton = false;
spectrogram_parameters.verbose = true;

spectrogram_parameters.weighting = 'unity';
[ spect, stimes, sfreqs ] = multitaper_spectrogram_mex(single(EEG), Fs, ...
    [0 min([Fs/2 spectrogram_parameters.frequency_max])], ...
    spectrogram_parameters.taper_params, spectrogram_parameters.window_params, ...
    spectrogram_parameters.min_NFFT, spectrogram_parameters.detrend, spectrogram_parameters.weighting,...
    spectrogram_parameters.ploton, spectrogram_parameters.verbose);


[ spect1, stimes, sfreqs ] = multitaper_spectrogram(single(EEG), Fs, ...
    [0 min([Fs/2 spectrogram_parameters.frequency_max])], ...
    spectrogram_parameters.taper_params, spectrogram_parameters.window_params, ...
    spectrogram_parameters.min_NFFT, spectrogram_parameters.detrend, spectrogram_parameters.weighting,...
    spectrogram_parameters.ploton, spectrogram_parameters.verbose);


% figure 
% ax1=subplot(2,1,1);
% imagesc(stimes, sfreqs, pow2db(spect))
% axis xy
% colormap jet 
% cax = climscale;
% title('Mex')
% 
% ax2=subplot(2,1,2);
% imagesc(stimes, sfreqs, pow2db(spect1))
% axis xy
% colormap jet 
% caxis(cax)
% title('Reg')
% 
% 
% linkaxes([ax1,ax2], 'xy')

mean(abs(spect(:) - spect1(:)))
std(abs(spect(:) - spect1(:)))

ratios = spect(:) ./ spect1(:);

scatter([0:1/Fs:(9999/Fs)],randsample(ratios,10000));
