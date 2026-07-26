import sys
from tkinter import Tk, Label, Button, GROOVE, CENTER, StringVar, END, LEFT, ttk, filedialog
from PIL import ImageTk, Image
import os
import pyedflib
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.signal.windows import dpss
from scipy.signal import detrend
from scipy.fft import rfft as scipy_rfft  # Tier D: multithreaded FFT
import warnings
import timeit
from functools import partial
from multiprocessing import Pool, cpu_count


# Helper #
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)


# Helper #
def center(toplevel):
    toplevel.update_idletasks()

    # Tkinter way to find the screen resolution
    screen_width = toplevel.winfo_screenwidth()
    screen_height = toplevel.winfo_screenheight()

    size = tuple(int(_) for _ in toplevel.geometry().split('+')[0].split('x'))
    x = screen_width/2 - size[0]/2
    y = screen_height/2 - size[1]/2

    toplevel.geometry("+%d+%d" % (x, y))


# Helper #
def change_entry_type(entry, entry_name, desired_type):
    if desired_type == "int":
        try:
            out = int(entry)
        except ValueError:
            print(entry_name + " could not be converted to an integer")
            out = 0
    elif desired_type == "float":
        try:
            out = float(entry)
        except ValueError:
            print(entry_name + " could not be converted to a float")
            out = 0
    else:
        raise ValueError('desired output not supported')

    return out


# Helper #
def is_outlier(data):
    smad = 1.4826 * np.median(abs(data - np.median(data)))  # scaled median absolute deviation
    outlier_mask = abs(data-np.median(data)) > 3*smad  # outliers are more than 3 smads away from median
    outlier_mask = (outlier_mask | np.isnan(data) | np.isinf(data))
    return outlier_mask


# Helper #
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

def plot_spect(root):
    '''
    '''

    # Plot Spectrogram with scaled colors #
    # Eliminate outliers and bad data from colormap scaling
    outlier_mask = np.apply_along_axis(is_outlier, 1, root.spect)  # apply outlier function to each column
    spect_data = root.spect[~outlier_mask]
    clim = np.percentile(spect_data, [5, 98])  # Scale colormap from 5th percentile to 98th

    # Plot #
    plt.figure(1, figsize=(10, 5))

    plt.pcolormesh(root.stimes, root.sfreqs, nanpow2db(root.spect), shading='auto', cmap='jet')
    plt.colorbar()
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.show()


def download_spect(root):
    '''
    '''
    root.saveloc = filedialog.askdirectory()
    np.savetxt(os.path.join(root.saveloc, 'spect.csv'), root.spect, delimiter=',')
    np.savetxt(os.path.join(root.saveloc, 'stimes.csv'), root.stimes, delimiter=',')
    np.savetxt(os.path.join(root.saveloc, 'sfreqs.csv'), root.sfreqs, delimiter=',')


def go_multitaper_spectrogram(root):
    """

    :param root:
    :return:
    """

    # Get all multitaper params and convert to appropriate data types #
    data = root.data
    fs = root.fs
    frequency_range = [change_entry_type(root.freq_from_entered.get(), "Frequency Min", "float"),
                       change_entry_type(root.freq_to_entered.get(), "Frequency Max", "float")]
    time_bandwidth = change_entry_type(root.hw_entered.get(), "Time Half Bandwidth", "int")
    num_tapers = change_entry_type(root.num_tapers_entered.get(), "Number of Tapers", "int")
    window_params = [change_entry_type(root.w_size_entered.get(), "Window Size", "int"),
                     change_entry_type(root.w_step_entered.get(), "Window Step Size", "int")]
    min_nfft = change_entry_type(root.min_nfft_entered.get(), "Minimum NFFT", "int")
    detrend = root.detrend_method_entered.get()
    multiprocess = True
    cpus = False
    plot_on = False
    weighting = root.weighting_entered.get()

    # Run multitaper spectrogram #
    root.spect, root.stimes, root.sfreqs = multitaper_spectrogram(data, fs, frequency_range, time_bandwidth, num_tapers,
                                                                  window_params, min_nfft, detrend, multiprocess, cpus,
                                                                  weighting, plot_on)

    # Make button for downloading spectrogram data
    download_data_button = Button(root, text="Download Output .csv", width=20,
                               height=2, relief=GROOVE, bg='#0064BB', fg="white", highlightbackground='#0064BB',
                               command=lambda: download_spect(root))
    download_data_button.place(relx=0.7, rely=0.92)

    plot_button = Button(root, text="Plot Spectrogram", width=20,
                               height=2, relief=GROOVE, bg='#0064BB', fg="white", highlightbackground='#0064BB',
                               command=lambda: plot_spect(root))
    plot_button.place(relx=0.15, rely=0.92)




def param_info_popup():
    param_info = Tk()
    param_info.geometry("500x650")
    param_info.title('Multitaper Parameter Info')
    descriptions = Label(param_info, text="*Frequency Range*: Range of frequencies (Hz) across which to compute the "
                                          "spectrum (default=0-50Hz).\n\n*Half-time Bandwidth (TW)*: TW can be computed as N*(BW/2) where "
                                          "N is the length of the window (seconds) and BW is the bandwidth of the main "
                                          "lobe. The bandwidth of the main lobe is also called the frequency resolution "
                                          "because it dictates the minimum difference in frequency that can be detected"
                                          "(default=5).\n\n*Number of Tapers (L)*: The number of DPSS tapers to be used "
                                          "to compute the spectrum. The optimal number of tapers is 2*(TW)-1. The "
                                          "default is 9 tapers.\n\n*Window Size (N) and Window Step Size*: These "
                                          "parameters dictate the temporal resolution of the analysis. The multitaper "
                                          "spectrum is computed for a single window of data, then the window moves "
                                          "based on step size and the spectrum will be computed again on the next window "
                                          "until the whole data array has been covered. Default window size is 10s, "
                                          "while default window step size is 5s.\n\n*Window Detrend Method*: Each window "
                                          "of data can be detrended to remove very low frequency oscillation artifacts "
                                          "that can come from a variety of sources. In 'linear' detrending, a linear "
                                          "model is fit to the window then subtracted out, while in 'constant' "
                                          "detrending the data is set to be zero mean. The default is 'linear', and the "
                                          "options are 'linear', 'constant', and 'off'.)\n\n*Minimum NFFT*: Multitaper "
                                          "spectrum computation relies on the Fourier Transform to transform the data "
                                          "from the time domain to the frequency domain. The Fast Fourier Transform "
                                          "(FFT) is an very computationally efficient algorithm to compute the Fourier "
                                          "Transform, and it works most efficiently when the number of data points in "
                                          "the given time series is a power of 2. Therefore, we want to pad the data "
                                          "with 0s to make it reach the closest power of 2. This implementation pads "
                                          "with zeros to the nearest power of 2 automatically, but if a specific power "
                                          "of 2 above the closest power fo 2 is desired, use this parameter. The default "
                                          "is 0.\n\n*Taper Weighting Method*: The DPSS tapers can be weighted "
                                          "differently, and we have included 2 weighting method options - 'eigen' and "
                                          "'adaptive' - along with the uniformly weighted option 'unity' which is the "
                                          "default for all implementations. Eigenvalue weighting weights the "
                                          "contribution of each taper to the spectrum by it's eigenvalue (frequency "
                                          "concentration). In most cases this makes little difference because most "
                                          "taper's eigenvalues are very close to one anyway. The adaptive weighting "
                                          "method weights the tapers in such a way as to reduce the broadband leakage "
                                          "of non-white ('colored') noise. This method is adapted from pages 368-370 of "
                                          "Percival and Walden's 'Spectral Analysis for Physical Applications: "
                                          "Multitaper and Conventional Univariate Techniques'. In practice, the adaptive "
                                          "method does not change the results much at all but is provided here for the "
                                          "sake of completeness.",
                         justify=LEFT, wraplength=500, padx=10, pady=10)
    descriptions.pack()


def channel_changed(event, root):
    """

    :param event:
    :param root:
    :return:
    """

    # Read edf file
    loading_label = Label(root, text="Loading in EDF channel...")
    loading_label.place(relx=0.66, rely=0.45, anchor=CENTER)
    root.update()
    root.signals, root.signal_headers, root.header = pyedflib.highlevel.read_edf(root.filename)
    loading_label.configure(text='Channel Loaded')

    # Get data and fs from extracted edf data #
    channel = root.channel.get()
    index_channel = [i for i in range(len(root.signal_headers)) if root.signal_headers[i]['label'] == channel][0]
    root.fs = root.signal_headers[index_channel]['sample_rate']
    root.data = root.signals[index_channel]
    Label(root, text="Sampling frequency of selected channel: "
                     + str(root.fs) + " Hz").place(relx=0.15, rely=0.55, anchor=CENTER)

    # Populate other parameters #
    # Freq range param
    y = 0.6
    Label(root, text="Frequency Range of Analysis:").place(relx=0.1, rely=y, anchor=CENTER)
    root.freq_from_entered = ttk.Entry(root, width=10)
    root.freq_from_entered.place(relx=0.225, rely=y, anchor=CENTER)
    root.freq_from_entered.delete(0, END)
    root.freq_from_entered.insert(0, 0)
    Label(root, text=" - ").place(relx=0.275, rely=y, anchor=CENTER)
    root.freq_to_entered = ttk.Entry(root, width=10)
    root.freq_to_entered.place(relx=0.325, rely=y, anchor=CENTER)
    root.freq_to_entered.delete(0, END)
    root.freq_to_entered.insert(0, 50)
    Label(root, text=" (Hz)").place(relx=0.375, rely=y, anchor=CENTER)

    # Taper Parameters
    y = 0.64
    Label(root, text="Half-Time Bandwidth (TW):").place(relx=0.096, rely=y, anchor=CENTER)
    root.hw_entered = ttk.Entry(root, width=10)
    root.hw_entered.delete(0, END)
    root.hw_entered.insert(0, 5)
    root.hw_entered.place(relx=0.22, rely=y, anchor=CENTER)

    y = 0.68
    Label(root, text="Number of Tapers (L):").place(relx=0.08, rely=y, anchor=CENTER)
    root.num_tapers_entered = ttk.Entry(root, width=10)
    root.num_tapers_entered.place(relx=0.19, rely=y, anchor=CENTER)
    root.num_tapers_entered.delete(0, END)
    root.num_tapers_entered.insert(0, 9)

    # Window Parameters
    y = 0.72
    Label(root, text="Window Size (N):").place(relx=0.0675, rely=y, anchor=CENTER)
    root.w_size_entered = ttk.Entry(root, width=10)
    root.w_size_entered.place(relx=0.16, rely=y, anchor=CENTER)
    root.w_size_entered.delete(0, END)
    root.w_size_entered.insert(0, 10)
    Label(root, text="(s)").place(relx=0.195, rely=y, anchor=CENTER)

    y = 0.76
    Label(root, text="Window Step Size:").place(relx=0.072, rely=y, anchor=CENTER)
    root.w_step_entered = ttk.Entry(root, width=10)
    root.w_step_entered.place(relx=0.16, rely=y, anchor=CENTER)
    root.w_step_entered.delete(0, END)
    root.w_step_entered.insert(0, 5)
    Label(root, text="(s)").place(relx=0.185, rely=y, anchor=CENTER)

    # Detrend method param
    y = 0.8
    Label(root, text="Window Detrend Method: ").place(relx=0.094, rely=y, anchor=CENTER)
    root.detrend_method_entered = ttk.Entry(root, width=10)
    root.detrend_method_entered.place(relx=0.2, rely=y, anchor=CENTER)
    root.detrend_method_entered.delete(0, END)
    root.detrend_method_entered.insert(0, 'linear')
    Label(root, text="Choices: 'linear' - default, 'constant', 'off'").place(relx=0.345, rely=y, anchor=CENTER)

    # Minimum NFFT Param
    y = 0.84
    Label(root, text="Minimum NFFT:").place(relx=0.069, rely=y, anchor=CENTER)
    root.min_nfft_entered = ttk.Entry(root, width=10)
    root.min_nfft_entered.place(relx=0.15, rely=y, anchor=CENTER)
    root.min_nfft_entered.delete(0, END)
    root.min_nfft_entered.insert(0, 0)

    # Taper Weighting Param
    y = 0.88
    Label(root, text="Taper Weighting Method:").place(relx=0.09, rely=y, anchor=CENTER)
    root.weighting_entered = ttk.Entry(root, width=10)
    root.weighting_entered.place(relx=0.2, rely=y, anchor=CENTER)
    root.weighting_entered.delete(0, END)
    root.weighting_entered.insert(0, 'unity')
    Label(root, text="Choices: 'unity' - default/recommended, 'eigen', 'adapt'").place(relx=0.385,
                                                                                       rely=y, anchor=CENTER)

    # Button for more parameter information
    param_info_button = Button(root, text="Parameter Descriptions", width=20,
                               height=2, relief=GROOVE, bg='#0064BB', fg="white", highlightbackground='#0064BB',
                               command=param_info_popup)
    param_info_button.place(relx=0.02, rely=0.46)

    # Multitaper Calc Button
    multitaper_go = Button(root, text="Calculate Multitaper Spectrogram", width=25, height=2, relief=GROOVE,
                           bg='#0064BB', fg="white", highlightbackground='#0064BB',
                           command=lambda: go_multitaper_spectrogram(root))
    multitaper_go.place(relx=0.5, rely=0.95, anchor=CENTER)


def start_edf(root):
    """

    :param root:
    :return:
    """

    # Load the file
    root.filename = filedialog.askopenfilenames(initialdir="~/", title="Select EDF file",
                                                filetypes=(("EDF files", "*.edf"), ("All files", "*.*")))[0]

    # Get channel names from edf
    print(root.filename)
    f = pyedflib.EdfReader(root.filename)
    root.signal_labels = f.getSignalLabels()

    # Create text label for select channel dropdown
    label_channel = Label(root, text='Select Channel:')
    label_channel.place(relx=0.39, rely=0.45, anchor=CENTER)

    # Create select channel dropdown
    selected_channel = StringVar()
    root.channel = ttk.Combobox(root, textvariable=selected_channel)
    root.channel['values'] = root.signal_labels
    root.channel['state'] = 'readonly'
    root.channel.place(relx=0.52, rely=0.45, anchor=CENTER)
    root.channel.bind('<<ComboboxSelected>>', lambda event: channel_changed(event, root))

    return root.filename


def main():
    # Create window
    root = Tk()
    root.geometry("1000x700")
    root.title('Multitaper Spectrogram Analysis')
    center(root)

    # Put image in window
    spath = resource_path("multitaper_image.png")
    simg = Image.open(spath)
    simg = simg.resize((700, 200))
    simg = ImageTk.PhotoImage(simg)
    panel = Label(root, image=simg)
    center(root)
    panel.pack()

    # Put button on window
    btn = Button(root, text="Choose EDF File", width=25, height=2, relief=GROOVE, bg='#0064BB',
                 fg="white", command=lambda: start_edf(root), highlightbackground='#0064BB')
    btn.place(relx=0.5, rely=0.39, anchor=CENTER)

    choose_file = Label(root, text="Select EDF file to begin")
    choose_file.place(relx=0.5, rely=0.34, anchor=CENTER)

    root.mainloop()


# MULTITAPER SPECTROGRAM #
def multitaper_spectrogram(data, fs, frequency_range=None, time_bandwidth=5, num_tapers=None, window_params=None,
                           min_nfft=0, detrend_opt='linear', multiprocess=False, cpus=False, weighting='unity',
                           plot_on=True, verbose=True, xyflip=False):
    """ Compute multitaper spectrogram of timeseries data
    Usage:
    mt_spectrogram, stimes, sfreqs = multitaper_spectrogram(data, fs, frequency_range=None, time_bandwidth=5,
                                                                   num_tapers=None, window_params=None, min_nfft=0,
                                                                   detrend_opt='linear', multiprocess=False, cpus=False,
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
                cpus (int): Number of cpus to use if multiprocess = True (default: False). Note: if default is left
                            as False and multiprocess = True, the number of cpus used for multiprocessing will be
                            all available - 1.
                weighting (str): weighting of tapers ('unity' (default), 'eigen', 'adapt');
                plot_on (bool): plot results (default: True)
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
            cpus = 3  # use 3 cores in multiprocessing
            weighting = 'unity'  # weight each taper at 1
            plot_on = True  # plot spectrogram
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
                                                           cpus, weighting, plot_on, verbose, xyflip):

        This code is companion to the paper:
        "Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"
           Michael J. Prerau, Ritchie E. Brown, Matt T. Bianchi, Jeffrey M. Ellenbogen, Patrick L. Purdon
           December 7, 2016 : 60-92
           DOI: 10.1152/physiol.00062.2015
         which should be cited for academic use of this code.

         A full tutorial on the multitaper spectrogram can be found at:  #   http://www.sleepEEG.org/multitaper

        Copyright 2021 Michael J. Prerau Laboratory. - http://www.sleepEEG.org
        Authors: Michael J. Prerau, Ph.D., Thomas Possidente

        Last modified - 2/18/2021 Thomas Possidente
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

    # set all but 1 arg of calc_mts_segment to constant (so we only have to supply one argument later)
    calc_mts_segment_plus_args = partial(calc_mts_segment, dpss_tapers=dpss_tapers, dpss_tapers_T=dpss_tapers_T,
                                         nfft=nfft, freq_inds=freq_inds,
                                         detrend_opt=detrend_opt, num_tapers=num_tapers, dpss_eigen=dpss_eigen,
                                         weighting=weighting, wt=wt)

    # Tier D: batched FFT path for 'unity'/'eigen'; 'adapt' keeps per-window loop.
    # Batched path uses scipy.fft workers=-1 (threading), so `multiprocess` flag is
    # ignored for non-adapt weighting.
    if weighting == 'adapt':
        if multiprocess:
            if not cpus:
                pool = Pool(cpu_count() - 1)
            else:
                pool = Pool(cpus)
            mt_spectrogram = pool.map(calc_mts_segment_plus_args, data_segments)
            pool.close()
            pool.join()
        else:
            mt_spectrogram = np.apply_along_axis(calc_mts_segment_plus_args, 1, data_segments)
    else:
        mt_spectrogram = calc_mts_batch(data_segments, dpss_tapers_T, nfft, freq_inds,
                                        detrend_opt, weighting, wt)

    # Compute one-sided PSD spectrum
    mt_spectrogram = np.asarray(mt_spectrogram)
    mt_spectrogram = mt_spectrogram.T
    dc_select = np.where(np.isclose(sfreqs, 0))
    nyquist_select = np.where(np.isclose(sfreqs, fs / 2))
    select = np.setdiff1d(np.arange(0, len(sfreqs)), [dc_select, nyquist_select])

    mt_spectrogram = np.vstack([mt_spectrogram[dc_select[0], :], 2 * mt_spectrogram[select, :],
                                mt_spectrogram[nyquist_select[0], :]]) / fs

    # Flip if requested
    if xyflip:
        mt_spectrogram = np.transpose(mt_spectrogram)

    # End timer and get elapsed compute time
    toc = timeit.default_timer()
    elapsed_time = toc - tic
    if verbose:
        print("\n Multitaper compute time: " + str(elapsed_time) + " seconds")

    # Plot multitaper spectrogram
    if plot_on:
        plt.figure(1, figsize=(10, 5))
        plt.pcolormesh(stimes, sfreqs, nanpow2db(mt_spectrogram), shading='auto', cmap='jet')
        plt.colorbar(label='Power (dB)')
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency (Hz)")
        plt.show()

    # Put outputs into better format for output
    # stimes = np.mat(stimes)
    # sfreqs = np.mat(sfreqs)

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
                    winsize_samples (float): number of samples in single time window
                    winstep_samples (float): number of samples in a single window step
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
        elif detrend_opt in ['none', 'false']:
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
             winsize_samples, winstep_samples, window_start, num_windows, nfft,
             detrend_opt, plot_on, verbose])


# PROCESS THE SPECTROGRAM PARAMETERS #
def process_spectrogram_params(fs, nfft, frequency_range, window_start, datawin_size):
    """ Helper function to create frequency vector and window indices
        Arguments:
             fs (float): sampling frequency in Hz  -- required
             nfft (int): length of signal to calculate fft on -- required
             frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] -- required
             window_start (1xm np.array): array of timestamps representing the beginning time for each
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

    # If segment has all zeros, return vector of zeros (Tier 1b)
    if all(data_segment == 0):
        return np.zeros(int(freq_inds.sum()))

    # Option to detrend data to remove low frequency DC component
    if detrend_opt != 'off':
        data_segment = detrend(data_segment, type=detrend_opt)

    # Multiply data by dpss tapers (STEP 2) - Tier 1a: use precomputed transpose
    tapered_data = dpss_tapers_T * data_segment[:, np.newaxis]

    # Compute the FFT (STEP 3) - Tier 3a: rfft
    fft_data = np.fft.rfft(tapered_data, nfft, axis=0)

    # Compute the weighted mean spectral power across tapers (STEP 4) - Tier 2a
    r = fft_data.real
    i = fft_data.imag
    spower = r * r + i * i
    nfreq = spower.shape[0]  # Tier 3f: one-sided length
    if weighting == 'adapt':
        # adaptive weights - for colored noise spectrum (Percival & Walden p368-370)
        tpower = np.dot(np.transpose(data_segment), (data_segment / len(data_segment)))
        spower_iter = np.mean(spower[:, 0:2], 1)
        spower_iter = spower_iter[:, np.newaxis]
        a = (1 - dpss_eigen) * tpower
        for i in range(3):  # 3 iterations only
            # Calc the MSE weights
            b = np.dot(spower_iter, np.ones((1, num_tapers))) / ((np.dot(spower_iter, np.transpose(dpss_eigen))) +
                                                                 (np.ones((nfreq, 1)) * np.transpose(a)))
            # Calc new spectral estimate
            wk = (b ** 2) * np.dot(np.ones((nfreq, 1)), np.transpose(dpss_eigen))
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
    """Vectorized multitaper spectrum over many windows (chunked batched FFT).

    See main-module docstring for full description. Used for 'unity' and 'eigen'
    weighting; 'adapt' keeps the per-window path.
    """
    W = data_segments.shape[0]
    nfreq_out = int(freq_inds.sum())
    out = np.empty((W, nfreq_out), dtype=np.float64)

    for start in range(0, W, batch_size):
        end = min(start + batch_size, W)
        batch = data_segments[start:end]
        chunk_out = out[start:end]

        # Match the main module and MATLAB behavior: any NaN invalidates the
        # whole window, which must not be passed to scipy.signal.detrend.
        valid_mask = ~np.any(np.isnan(batch), axis=1)
        if not np.all(valid_mask):
            chunk_out.fill(np.nan)
            if not np.any(valid_mask):
                continue
            batch = batch[valid_mask]

        zero_mask = ~np.any(batch, axis=1)

        if detrend_opt != 'off':
            batch = detrend(batch, type=detrend_opt, axis=1)

        tapered = batch[:, :, None] * dpss_tapers_T[None, :, :]
        fft_data = scipy_rfft(tapered, n=nfft, axis=1, workers=-1)
        spower = fft_data.real ** 2 + fft_data.imag ** 2

        if weighting == 'unity':
            mt = spower.mean(axis=2)
        else:
            mt = np.squeeze(spower @ wt, axis=-1)

        mt_filtered = mt[:, freq_inds]
        if np.any(zero_mask):
            mt_filtered[zero_mask, :] = 0

        chunk_out[valid_mask] = mt_filtered

    return out


if __name__ == '__main__':
    main()




# TODO: loading bar for ingesting edfs, loading icon/bar for spectrogram calc
