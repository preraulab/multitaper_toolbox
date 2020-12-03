function [mt_spectrogram,stimes,sfreqs] = multitaper_spectrogram_optimized(varargin)
    [mt_spectrogram,stimes,sfreqs] = multitaper_spectrogram(varargin{:});
    warning('multitaper_spectrogram_optimized.m is deprecated. Please call multitaper_spectrogram.m instead.')
end
