% Compare multitaper estimation methods 
close all; clear all;

% set up simualtion data 
duration = 60; %60s
Fs = 500; %100Hz
EEG = single(pinknoise(Fs*duration))';
nw = 2;

% MTM method 
spectrogram_parameters.frequency_max = Fs/2;
spectrogram_parameters.taper_params = [nw, 2*nw-1];
spectrogram_parameters.window_params = [1, 0.05];
spectrogram_parameters.min_NFFT = 2^10;
spectrogram_parameters.detrend = 'constant';
spectrogram_parameters.weighting = 'unity';
spectrogram_parameters.ploton = false;
spectrogram_parameters.verbose = true;

[ spect, stimes, sfreqs ] = multitaper_spectrogram_mex(EEG, Fs, ...
    [0 min([Fs/2 spectrogram_parameters.frequency_max])], ...
    spectrogram_parameters.taper_params, spectrogram_parameters.window_params, ...
    spectrogram_parameters.min_NFFT, spectrogram_parameters.detrend, spectrogram_parameters.weighting,...
    spectrogram_parameters.ploton, spectrogram_parameters.verbose);

spectrogram_parameters.weighting = 'eigen';
[ spect1, stimes, sfreqs ] = multitaper_spectrogram_mex(EEG, Fs, ...
    [0 min([Fs/2 spectrogram_parameters.frequency_max])], ...
    spectrogram_parameters.taper_params, spectrogram_parameters.window_params, ...
    spectrogram_parameters.min_NFFT, spectrogram_parameters.detrend, spectrogram_parameters.weighting,...
    spectrogram_parameters.ploton, spectrogram_parameters.verbose);

spectrogram_parameters.weighting = 'adapt';
[ spect2, stimes, sfreqs ] = multitaper_spectrogram_mex(EEG, Fs, ...
    [0 min([Fs/2 spectrogram_parameters.frequency_max])], ...
    spectrogram_parameters.taper_params, spectrogram_parameters.window_params, ...
    spectrogram_parameters.min_NFFT, spectrogram_parameters.detrend, spectrogram_parameters.weighting,...
    spectrogram_parameters.ploton, spectrogram_parameters.verbose);

figure 
ax1=subplot(3,1,1);
imagesc(stimes, sfreqs, pow2db(spect))
axis xy
colormap jet 
cax = climscale;
title('Unity')

ax2=subplot(3,1,2);
imagesc(stimes, sfreqs, pow2db(spect1))
axis xy
colormap jet 
caxis(cax)
title('Eigen')

ax3=subplot(3,1,3);
imagesc(stimes, sfreqs, pow2db(spect2))
axis xy
colormap jet 
caxis(cax)
title('Adapt')

linkaxes([ax1,ax2,ax3], 'xy')


% MATLAB pmtm method 
[ pxx, f ] = pmtm(EEG, nw, [], Fs, 'adapt');

% compare the two methods 
figure
plot(sfreqs, pow2db(spect), 'LineWidth', 1)
hold on
plot(f, pow2db(pxx), 'LineWidth', 1)
legend('Multitaper', 'pmtm')

figure
plot(pxx./spect)