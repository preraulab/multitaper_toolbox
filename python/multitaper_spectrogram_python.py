# Analysis Imports
import math
import sys
import numpy as np
from scipy.signal.windows import dpss
from scipy.signal import detrend
from scipy.fft import rfft as scipy_rfft  # Tier D: multithreaded FFT
# Logistical Imports
import warnings
import timeit
from joblib import Parallel, delayed, cpu_count
# Visualization imports
import matplotlib.pyplot as plt

# Optional Rust backend. If the compiled extension is unavailable we fall
# back transparently to the pure-Python path. This import must never fail
# the module load on machines without the Rust extension installed.
try:
    from multitaper_rs import compute_spectrogram as _rust_compute
    _HAS_RUST = True
except ImportError:
    _rust_compute = None
    _HAS_RUST = False

# Guard so the "which backend" banner only prints once per process.
_BACKEND_BANNER_PRINTED = False


def _announce_backend(using_rust):
    """Print a one-time banner to stderr indicating which backend is in use."""
    global _BACKEND_BANNER_PRINTED
    if _BACKEND_BANNER_PRINTED:
        return
    _BACKEND_BANNER_PRINTED = True
    if using_rust:
        print("multitaper: using Rust backend", file=sys.stderr)
    elif _HAS_RUST:
        print("multitaper: Rust backend available but disabled, using pure Python",
              file=sys.stderr)
    else:
        print("multitaper: Rust backend unavailable, using pure Python",
              file=sys.stderr)


# MULTITAPER SPECTROGRAM #
def multitaper_spectrogram(data, fs, frequency_range=None, time_bandwidth=5, num_tapers=None, window_params=None,
                           min_nfft=0, detrend_opt='linear', multiprocess=False, n_jobs=None, weighting='unity',
                           plot_on=True, clim_scale=True, verbose=True, xyflip=False, use_rust=None):
    """ Compute multitaper spectrogram of timeseries data
    Usage:
    mt_spectrogram, stimes, sfreqs = multitaper_spectrogram(data, fs, frequency_range=None, time_bandwidth=5,
                                                                   num_tapers=None, window_params=None, min_nfft=0,
                                                                   detrend_opt='linear', multiprocess=False, n_jobs=None,
                                                                    weighting='unity', plot_on=True, verbose=True,
                                                                    xyflip=False):
        Arguments:
                data (1d np.array): time series data -- required
                fs (float): sampling frequency in Hz  -- required
                frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] (default: [0 nyquist])
                time_bandwidth (float): time-half bandwidth product (window duration*half bandwidth of main lobe)
                                        (default: 5 Hz*s)
                num_tapers (int): number of DPSS tapers to use (default: [will be computed
                                  as floor(2*time_bandwidth - 1)])
                window_params (list): 1x2 list - [window size (seconds), step size (seconds)] (default: [5 1])
                detrend_opt (string): detrend data window ('linear' (default), 'constant', 'off')
                                      (Default: 'linear')
                min_nfft (int): minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x)
                                (default: 0)
                multiprocess (bool): Use multiprocessing to compute multitaper spectrogram (default: False)
                n_jobs (int): Number of cpus to use if multiprocess = True (default: False). Note: if default is left
                            as None and multiprocess = True, the number of cpus used for multiprocessing will be
                            all available - 1.
                weighting (str): weighting of tapers ('unity' (default), 'eigen', 'adapt');
                plot_on (bool): plot results (default: True)
                clim_scale (bool): automatically scale the colormap on the plotted spectrogram (default: true)
                verbose (bool): display spectrogram properties (default: true)
                xyflip (bool): transpose the mt_spectrogram output (default: false)
        Returns:
                mt_spectrogram (TxF np array): spectral power matrix
                stimes (1xT np array): timepoints (s) in mt_spectrogram
                sfreqs (1xF np array)L frequency values (Hz) in mt_spectrogram

        Example:
        In this example we create some chirp data and run the multitaper spectrogram on it.
            import numpy as np  # import numpy
            from scipy.signal import chirp  # import chirp generation function
            # Set spectrogram params
            fs = 200  # Sampling Frequency
            frequency_range = [0, 25]  # Limit frequencies from 0 to 25 Hz
            time_bandwidth = 3  # Set time-half bandwidth
            num_tapers = 5  # Set number of tapers (optimal is time_bandwidth*2 - 1)
            window_params = [4, 1]  # Window size is 4s with step size of 1s
            min_nfft = 0  # No minimum nfft
            detrend_opt = 'constant'  # detrend each window by subtracting the average
            multiprocess = True  # use multiprocessing
            n_jobs = 3  # use 3 cores in multiprocessing
            weighting = 'unity'  # weight each taper at 1
            plot_on = True  # plot spectrogram
            clim_scale = False # don't auto-scale the colormap
            verbose = True  # print extra info
            xyflip = False  # do not transpose spect output matrix
            # Generate sample chirp data
            t = np.arange(1/fs, 600, 1/fs)  # Create 10 min time array from 1/fs to 600 stepping by 1/fs
            f_start = 1  # Set chirp freq range min (Hz)
            f_end = 20  # Set chirp freq range max (Hz)
            data = chirp(t, f_start, t[-1], f_end, 'logarithmic')
            # Compute the multitaper spectrogram
            spect, stimes, sfreqs = multitaper_spectrogram(data, fs, frequency_range, time_bandwidth, num_tapers,
                                                           window_params, min_nfft, detrend_opt, multiprocess,
                                                           n_jobs, weighting, plot_on, verbose, xyflip):

        This code is companion to the paper:
        "Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"
           Michael J. Prerau, Ritchie E. Brown, Matt T. Bianchi, Jeffrey M. Ellenbogen, Patrick L. Purdon
           December 7, 2016 : 60-92
           DOI: 10.1152/physiol.00062.2015
         which should be cited for academic use of this code.

         A full tutorial on the multitaper spectrogram can be found at:  #   https://www.sleepEEG.org/multitaper

        Copyright 2021 Michael J. Prerau Laboratory. - https://www.sleepEEG.org
        Authors: Michael J. Prerau, Ph.D., Thomas Possidente, Mingjian He

  __________________________________________________________________________________________________________________
    """

    #  Process user input
    [data, fs, frequency_range, time_bandwidth, num_tapers,
     winsize_samples, winstep_samples, window_start,
     num_windows, nfft, detrend_opt, plot_on, verbose] = process_input(data, fs, frequency_range, time_bandwidth,
                                                                       num_tapers, window_params, min_nfft,
                                                                       detrend_opt, plot_on, verbose)

    # Set up spectrogram parameters
    [window_idxs, stimes, sfreqs, freq_inds] = process_spectrogram_params(fs, nfft, frequency_range, window_start,
                                                                          winsize_samples)
    # Display spectrogram parameters
    if verbose:
        display_spectrogram_props(fs, time_bandwidth, num_tapers, [winsize_samples, winstep_samples], frequency_range,
                                  detrend_opt)

    # Split data into segments and preallocate
    data_segments = data[window_idxs]

    # COMPUTE THE MULTITAPER SPECTROGRAM
    #     STEP 1: Compute DPSS tapers based on desired spectral properties
    #     STEP 2: Multiply the data segment by the DPSS Tapers
    #     STEP 3: Compute the spectrum for each tapered segment
    #     STEP 4: Take the mean of the tapered spectra

    # Compute DPSS tapers (STEP 1)
    dpss_tapers, dpss_eigen = dpss(winsize_samples, time_bandwidth, num_tapers, return_ratios=True)
    dpss_eigen = np.reshape(dpss_eigen, (num_tapers, 1))

    # pre-compute weights
    if weighting == 'eigen':
        wt = dpss_eigen / num_tapers
    elif weighting == 'unity':
        wt = np.ones(num_tapers) / num_tapers
        wt = np.reshape(wt, (num_tapers, 1))  # reshape as column vector
    else:
        wt = 0

    tic = timeit.default_timer()  # start timer

    # Precompute transpose of tapers (used every segment) - Tier 1a
    dpss_tapers_T = dpss_tapers.T

    # Decide whether to use the Rust backend. 'adapt' weighting is not yet
    # implemented in Rust so it always falls back. use_rust=False forces the
    # pure-Python path (for A/B testing / equivalence checks).
    if use_rust is None:
        use_rust_resolved = _HAS_RUST and weighting in ('unity', 'eigen')
    elif use_rust:
        if not _HAS_RUST:
            raise RuntimeError("use_rust=True was requested but the multitaper_rs "
                               "extension is not installed. Build it with "
                               "`cd multitaper/rust && maturin develop --release`.")
        if weighting == 'adapt':
            raise NotImplementedError("Rust backend does not support weighting='adapt'. "
                                      "Use use_rust=False or weighting='unity'/'eigen'.")
        use_rust_resolved = True
    else:
        use_rust_resolved = False

    _announce_backend(use_rust_resolved)

    if use_rust_resolved:
        # Rust path: pass the raw signal + tapers; Rust does the windowing,
        # detrend, FFT, taper-weighted power sum, and one-sided PSD scaling.
        # (DPSS itself is still computed by scipy above because porting the
        # Slepian eigensolver to Rust would require pulling in LAPACK.)
        winsize_samples_py = int(window_params[0] * fs) if window_params is not None else int(5 * fs)
        # Recover actual window_params in seconds that Rust expects. We trust
        # the (possibly adjusted) winsize_samples from process_input.
        winsize_s = dpss_tapers.shape[1] / fs
        # winstep in seconds = winstep_samples / fs; recover from window_start spacing.
        if len(window_start) >= 2:
            winstep_s = (window_start[1] - window_start[0]) / fs
        elif window_params is not None:
            winstep_s = window_params[1]
        else:
            winstep_s = 1.0

        eigen_arr = None
        if weighting == 'eigen':
            eigen_arr = np.ascontiguousarray(dpss_eigen.reshape(-1), dtype=np.float64)

        data_c = np.ascontiguousarray(data, dtype=np.float64)
        tapers_c = np.ascontiguousarray(dpss_tapers, dtype=np.float64)

        mt_spectrogram, stimes_rs, sfreqs_rs = _rust_compute(
            data_c,
            tapers_c,
            float(fs),
            (float(frequency_range[0]), float(frequency_range[1])),
            (float(winsize_s), float(winstep_s)),
            int(nfft),
            detrend_opt,
            weighting,
            eigen_arr,
        )
        # Rust already applies one-sided PSD scaling and returns
        # (nfreq_out, num_windows). Replace the Python-computed stimes/sfreqs
        # with Rust's (they are algorithmically identical but keep alignment).
        stimes = stimes_rs
        sfreqs = sfreqs_rs

        if xyflip:
            mt_spectrogram = mt_spectrogram.T

        toc = timeit.default_timer()
        if verbose:
            print("\n Multitaper compute time (Rust): " + "%.2f" % (toc - tic) + " seconds")

        if plot_on:
            spect_data = mt_spectrogram
            clim = np.percentile(spect_data, [5, 95])
            plt.figure(1, figsize=(10, 5))
            plt.pcolormesh(stimes, sfreqs, nanpow2db(mt_spectrogram), shading='auto', cmap='jet')
            plt.colorbar(label='Power (dB)')
            plt.xlabel("Time (s)")
            plt.ylabel("Frequency (Hz)")
            if clim_scale:
                plt.clim(clim)
            plt.show()

        if np.all(mt_spectrogram.flatten() == 0):
            print("\n Data was all zeros, no output")

        return mt_spectrogram, stimes, sfreqs

    # Set up calc_mts_segment() input arguments
    mts_params = (dpss_tapers, dpss_tapers_T, nfft, freq_inds, detrend_opt, num_tapers, dpss_eigen, weighting, wt)

    # Tier D: batched FFT path for 'unity'/'eigen'; 'adapt' keeps the per-window loop
    # because its iterative convergence is per-window. Batched path uses scipy.fft
    # with workers=-1 (threading), which scales without joblib's process overhead,
    # so the `multiprocess` flag is ignored for non-adapt weighting.
    if weighting == 'adapt':
        if multiprocess:
            n_jobs = max(cpu_count() - 1, 1) if n_jobs is None else n_jobs
            mt_spectrogram = np.vstack(Parallel(n_jobs=n_jobs)(delayed(calc_mts_segment)(
                data_segments[num_window, :], *mts_params) for num_window in range(num_windows)))
        else:
            mt_spectrogram = np.apply_along_axis(calc_mts_segment, 1, data_segments, *mts_params)
    else:
        mt_spectrogram = calc_mts_batch(data_segments, dpss_tapers_T, nfft, freq_inds,
                                        detrend_opt, weighting, wt)

    # Compute one-sided PSD spectrum
    mt_spectrogram = mt_spectrogram.T
    dc_select = np.where(np.isclose(sfreqs, 0))[0]
    nyquist_select = np.where(np.isclose(sfreqs, fs/2))[0]
    select = np.setdiff1d(np.arange(0, len(sfreqs)), np.concatenate((dc_select, nyquist_select)))

    mt_spectrogram = np.vstack([mt_spectrogram[dc_select, :], 2*mt_spectrogram[select, :],
                               mt_spectrogram[nyquist_select, :]]) / fs

    # Flip if requested
    if xyflip:
        mt_spectrogram = mt_spectrogram.T

    # End timer and get elapsed compute time
    toc = timeit.default_timer()
    if verbose:
        print("\n Multitaper compute time: " + "%.2f" % (toc - tic) + " seconds")

    # Plot multitaper spectrogram
    if plot_on:

        # Eliminate bad data from colormap scaling
        spect_data = mt_spectrogram
        clim = np.percentile(spect_data, [5, 95])  # Scale colormap from 5th percentile to 95th

        plt.figure(1, figsize=(10, 5))
        plt.pcolormesh(stimes, sfreqs, nanpow2db(mt_spectrogram), shading='auto', cmap='jet')
        plt.colorbar(label='Power (dB)')
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency (Hz)")
        if clim_scale:
            plt.clim(clim)  # actually change colorbar scale
        plt.show()

    if all(mt_spectrogram.flatten() == 0):
        print("\n Data was all zeros, no output")

    return mt_spectrogram, stimes, sfreqs


# Helper Functions #

# Process User Inputs #
def process_input(data, fs, frequency_range=None, time_bandwidth=5, num_tapers=None, window_params=None, min_nfft=0,
                  detrend_opt='linear', plot_on=True, verbose=True):
    """ Helper function to process multitaper_spectrogram() arguments
            Arguments:
                    data (1d np.array): time series data-- required
                    fs (float): sampling frequency in Hz  -- required
                    frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] (default: [0 nyquist])
                    time_bandwidth (float): time-half bandwidth product (window duration*half bandwidth of main lobe)
                                            (default: 5 Hz*s)
                    num_tapers (int): number of DPSS tapers to use (default: None [will be computed
                                      as floor(2*time_bandwidth - 1)])
                    window_params (list): 1x2 list - [window size (seconds), step size (seconds)] (default: [5 1])
                    min_nfft (int): minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x)
                                    (default: 0)
                    detrend_opt (string): detrend data window ('linear' (default), 'constant', 'off')
                                          (Default: 'linear')
                    plot_on (True): plot results (default: True)
                    verbose (True): display spectrogram properties (default: true)
            Returns:
                    data (1d np.array): same as input
                    fs (float): same as input
                    frequency_range (list): same as input or calculated from fs if not given
                    time_bandwidth (float): same as input or default if not given
                    num_tapers (int): same as input or calculated from time_bandwidth if not given
                    winsize_samples (int): number of samples in single time window
                    winstep_samples (int): number of samples in a single window step
                    window_start (1xm np.array): array of timestamps representing the beginning time for each window
                    num_windows (int): number of windows in the data
                    nfft (int): length of signal to calculate fft on
                    detrend_opt ('string'): same as input or default if not given
                    plot_on (bool): same as input
                    verbose (bool): same as input
    """

    # Make sure data is 1 dimensional np array
    if len(data.shape) != 1:
        if (len(data.shape) == 2) & (data.shape[1] == 1):  # if it's 2d, but can be transferred to 1d, do so
            data = np.ravel(data[:, 0])
        elif (len(data.shape) == 2) & (data.shape[0] == 1):  # if it's 2d, but can be transferred to 1d, do so
            data = np.ravel(data.T[:, 0])
        else:
            raise TypeError("Input data is the incorrect dimensions. Should be a 1d array with shape (n,) where n is \
                            the number of data points. Instead data shape was " + str(data.shape))

    # Set frequency range if not provided
    if frequency_range is None:
        frequency_range = [0, fs / 2]

    # Set detrending method
    detrend_opt = detrend_opt.lower()
    if detrend_opt != 'linear':
        if detrend_opt in ['const', 'constant']:
            detrend_opt = 'constant'
        elif detrend_opt in ['none', 'false', 'off']:
            detrend_opt = 'off'
        else:
            raise ValueError("'" + str(detrend_opt) + "' is not a valid argument for detrend_opt. The choices " +
                             "are: 'constant', 'linear', or 'off'.")
    # Check if frequency range is valid
    if frequency_range[1] > fs / 2:
        frequency_range = list(frequency_range)
        frequency_range[1] = fs / 2
        warnings.warn('Upper frequency range greater than Nyquist, setting range to [' +
                      str(frequency_range[0]) + ', ' + str(frequency_range[1]) + ']')

    # Set number of tapers if none provided
    if num_tapers is None:
        num_tapers = math.floor(2 * time_bandwidth) - 1

    # Warn if number of tapers is suboptimal
    if num_tapers != math.floor(2 * time_bandwidth) - 1:
        warnings.warn('Number of tapers is optimal at floor(2*TW) - 1. consider using ' +
                      str(math.floor(2 * time_bandwidth) - 1))

    # If no window params provided, set to defaults
    if window_params is None:
        window_params = [5, 1]

    # Check if window size is valid, fix if not
    if window_params[0] * fs % 1 != 0:
        winsize_samples = round(window_params[0] * fs)
        warnings.warn('Window size is not divisible by sampling frequency. Adjusting window size to ' +
                      str(winsize_samples / fs) + ' seconds')
    else:
        winsize_samples = window_params[0] * fs

    # Check if window step is valid, fix if not
    if window_params[1] * fs % 1 != 0:
        winstep_samples = round(window_params[1] * fs)
        warnings.warn('Window step size is not divisible by sampling frequency. Adjusting window step size to ' +
                      str(winstep_samples / fs) + ' seconds')
    else:
        winstep_samples = window_params[1] * fs

    # Get total data length
    len_data = len(data)

    # Check if length of data is smaller than window (bad)
    if len_data < winsize_samples:
        raise ValueError("\nData length (" + str(len_data) + ") is shorter than window size (" +
                         str(winsize_samples) + "). Either increase data length or decrease window size.")

    # Find window start indices and num of windows
    window_start = np.arange(0, len_data - winsize_samples + 1, winstep_samples)
    num_windows = len(window_start)

    # Get num points in FFT
    if min_nfft == 0:  # avoid divide by zero error in np.log2(0)
        nfft = max(2 ** math.ceil(np.log2(abs(winsize_samples))), winsize_samples)
    else:
        nfft = max(max(2 ** math.ceil(np.log2(abs(winsize_samples))), winsize_samples),
                   2 ** math.ceil(np.log2(abs(min_nfft))))

    return ([data, fs, frequency_range, time_bandwidth, num_tapers,
             int(winsize_samples), int(winstep_samples), window_start, num_windows, nfft,
             detrend_opt, plot_on, verbose])


# PROCESS THE SPECTROGRAM PARAMETERS #
def process_spectrogram_params(fs, nfft, frequency_range, window_start, datawin_size):
    """ Helper function to create frequency vector and window indices
        Arguments:
             fs (float): sampling frequency in Hz  -- required
             nfft (int): length of signal to calculate fft on -- required
             frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] -- required
             window_start (1xm np array): array of timestamps representing the beginning time for each
                                          window -- required
             datawin_size (float): seconds in one window -- required
        Returns:
            window_idxs (nxm np array): indices of timestamps for each window
                                        (nxm where n=number of windows and m=datawin_size)
            stimes (1xt np array): array of times for the center of the spectral bins
            sfreqs (1xf np array): array of frequency bins for the spectrogram
            freq_inds (1d np array): boolean array of which frequencies are being analyzed in
                                      an array of frequencies from 0 to fs with steps of fs/nfft
    """

    # create frequency vector (Tier 3c: one-sided via rfftfreq)
    df = fs / nfft
    sfreqs = np.fft.rfftfreq(nfft, d=1/fs)

    # Get frequencies for given frequency range
    freq_inds = (sfreqs >= frequency_range[0]) & (sfreqs <= frequency_range[1])
    sfreqs = sfreqs[freq_inds]

    # Compute times in the middle of each spectrum
    window_middle_samples = window_start + round(datawin_size / 2)
    stimes = window_middle_samples / fs

    # Get indexes for each window
    window_idxs = np.atleast_2d(window_start).T + np.arange(0, datawin_size, 1)
    window_idxs = window_idxs.astype(int)

    return [window_idxs, stimes, sfreqs, freq_inds]


# DISPLAY SPECTROGRAM PROPERTIES
def display_spectrogram_props(fs, time_bandwidth, num_tapers, data_window_params, frequency_range, detrend_opt):
    """ Prints spectrogram properties
        Arguments:
            fs (float): sampling frequency in Hz  -- required
            time_bandwidth (float): time-half bandwidth product (window duration*1/2*frequency_resolution) -- required
            num_tapers (int): number of DPSS tapers to use -- required
            data_window_params (list): 1x2 list - [window length(s), window step size(s)] -- required
            frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] -- required
            detrend_opt (str): detrend data window ('linear' (default), 'constant', 'off')
        Returns:
            This function does not return anything
    """

    data_window_params = np.asarray(data_window_params) / fs

    # Print spectrogram properties
    print("Multitaper Spectrogram Properties: ")
    print('     Spectral Resolution: ' + str(2 * time_bandwidth / data_window_params[0]) + 'Hz')
    print('     Window Length: ' + str(data_window_params[0]) + 's')
    print('     Window Step: ' + str(data_window_params[1]) + 's')
    print('     Time Half-Bandwidth Product: ' + str(time_bandwidth))
    print('     Number of Tapers: ' + str(num_tapers))
    print('     Frequency Range: ' + str(frequency_range[0]) + "-" + str(frequency_range[1]) + 'Hz')
    print('     Detrend: ' + detrend_opt + '\n')


# NANPOW2DB
def nanpow2db(y):
    """ Power to dB conversion, setting bad values to nans
        Arguments:
            y (float or array-like): power
        Returns:
            ydB (float or np array): inputs converted to dB with 0s and negatives resulting in nans
    """

    if isinstance(y, int) or isinstance(y, float):
        if y == 0:
            return np.nan
        else:
            ydB = 10 * np.log10(y)
    else:
        if isinstance(y, list):  # if list, turn into array
            y = np.asarray(y)
        y = y.astype(float)  # make sure it's a float array so we can put nans in it
        y[y == 0] = np.nan
        ydB = 10 * np.log10(y)

    return ydB


# Helper #
def is_outlier(data):
    smad = 1.4826 * np.median(abs(data - np.median(data)))  # scaled median absolute deviation
    outlier_mask = abs(data-np.median(data)) > 3*smad  # outliers are more than 3 smads away from median
    outlier_mask = (outlier_mask | np.isnan(data) | np.isinf(data))
    return outlier_mask


# CALCULATE MULTITAPER SPECTRUM ON SINGLE SEGMENT
def calc_mts_segment(data_segment, dpss_tapers, dpss_tapers_T, nfft, freq_inds, detrend_opt, num_tapers, dpss_eigen, weighting, wt):
    """ Helper function to calculate the multitaper spectrum of a single segment of data
        Arguments:
            data_segment (1d np.array): One window worth of time-series data -- required
            dpss_tapers (2d np.array): Parameters for the DPSS tapers to be used.
                                       Dimensions are (num_tapers, winsize_samples) -- required
            nfft (int): length of signal to calculate fft on -- required
            freq_inds (1d np array): boolean array of which frequencies are being analyzed in
                                      an array of frequencies from 0 to fs with steps of fs/nfft
            detrend_opt (str): detrend data window ('linear' (default), 'constant', 'off')
            num_tapers (int): number of tapers being used
            dpss_eigen (np array):
            weighting (str):
            wt (int or np array):
        Returns:
            mt_spectrum (1d np.array): spectral power for single window
    """

    # If segment has all zeros, return vector of zeros (Tier 1b: use np.zeros)
    if all(data_segment == 0):
        return np.zeros(int(freq_inds.sum()))

    # Option to detrend data to remove low frequency DC component
    if detrend_opt != 'off':
        data_segment = detrend(data_segment, type=detrend_opt)

    # Multiply data by dpss tapers (STEP 2) - Tier 1a: use precomputed transpose
    tapered_data = dpss_tapers_T * data_segment[:, np.newaxis]

    # Compute the FFT (STEP 3) - Tier 3a: rfft (output length nfft//2+1)
    fft_data = np.fft.rfft(tapered_data, nfft, axis=0)

    # Compute the weighted mean spectral power across tapers (STEP 4) - Tier 2a
    r = fft_data.real
    i = fft_data.imag
    spower = r * r + i * i
    nfreq = spower.shape[0]  # Tier 3f: use spower length (one-sided)
    if weighting == 'adapt':
        # adaptive weights - for colored noise spectrum (Percival & Walden p368-370)
        tpower = np.dot(np.transpose(data_segment), (data_segment/len(data_segment)))
        spower_iter = np.mean(spower[:, 0:2], 1)
        spower_iter = spower_iter[:, np.newaxis]
        a = (1 - dpss_eigen) * tpower
        for i in range(3):  # 3 iterations only
            # Calc the MSE weights
            b = np.dot(spower_iter, np.ones((1, num_tapers))) / ((np.dot(spower_iter, np.transpose(dpss_eigen))) +
                                                                 (np.ones((nfreq, 1)) * np.transpose(a)))
            # Calc new spectral estimate
            wk = (b**2) * np.dot(np.ones((nfreq, 1)), np.transpose(dpss_eigen))
            spower_iter = np.sum((np.transpose(wk) * np.transpose(spower)), 0) / np.sum(wk, 1)
            spower_iter = spower_iter[:, np.newaxis]

        mt_spectrum = np.squeeze(spower_iter)

    elif weighting == 'unity':
        # Tier 2b: uniform weights == mean across tapers
        mt_spectrum = spower.mean(axis=1)
    else:
        # eigenvalue weights
        mt_spectrum = np.dot(spower, wt)
        mt_spectrum = np.reshape(mt_spectrum, nfreq)  # reshape to 1D

    return mt_spectrum[freq_inds]


# BATCHED MULTITAPER SEGMENT (Tier D) #
def calc_mts_batch(data_segments, dpss_tapers_T, nfft, freq_inds, detrend_opt, weighting, wt, batch_size=1024):
    """Vectorized multitaper spectrum over many windows at once.

    Replaces the per-window loop with chunked batched FFT (scipy.fft workers=-1).
    Used for 'unity' and 'eigen' weighting; 'adapt' keeps the per-window path
    because its iterative convergence is per-window.

    Arguments:
        data_segments (2d np.array): (W, winsize) — one row per window.
        dpss_tapers_T (2d np.array): (winsize, K) — precomputed taper transpose.
        nfft (int): FFT length (zero-padded).
        freq_inds (1d bool np.array): mask selecting in-range frequency bins
                                       (length nfft//2 + 1).
        detrend_opt (str): 'linear', 'constant', or 'off'.
        weighting (str): 'unity' or 'eigen'.
        wt (np.array): taper weights, shape (K, 1).
        batch_size (int): windows per chunk (bounds peak memory).

    Returns:
        mt_spectrogram (2d np.array): (W, freq_inds.sum()) — one row per window.
    """
    W = data_segments.shape[0]
    nfreq_out = int(freq_inds.sum())
    out = np.empty((W, nfreq_out), dtype=np.float64)

    for start in range(0, W, batch_size):
        end = min(start + batch_size, W)
        batch = data_segments[start:end]  # (chunk, winsize)

        # All-zero segment mask (captured before detrend, which cannot change all-zeros)
        zero_mask = ~np.any(batch, axis=1)

        # Row-wise detrend
        if detrend_opt != 'off':
            batch = detrend(batch, type=detrend_opt, axis=1)

        # Taper multiply: (chunk, winsize, 1) * (1, winsize, K) -> (chunk, winsize, K)
        tapered = batch[:, :, None] * dpss_tapers_T[None, :, :]

        # Batched multithreaded rfft along the time axis
        fft_data = scipy_rfft(tapered, n=nfft, axis=1, workers=-1)

        # Magnitude squared (Tier 2a form)
        spower = fft_data.real ** 2 + fft_data.imag ** 2  # (chunk, nfreq_full, K)

        # Taper-weighted mean
        if weighting == 'unity':
            mt = spower.mean(axis=2)
        else:  # 'eigen'
            mt = np.squeeze(spower @ wt, axis=-1)

        # Apply frequency-range mask and zero-segment override
        mt_filtered = mt[:, freq_inds]
        if np.any(zero_mask):
            mt_filtered[zero_mask, :] = 0

        out[start:end] = mt_filtered

    return out
