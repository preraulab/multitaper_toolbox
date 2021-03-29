from tkinter import *
from tkinter import ttk
from tkinter import filedialog

from PIL import ImageTk, Image

import sys
import os

import pyedflib
# import mne

import numpy as np

import matplotlib.pyplot as plt
import librosa.display

sys.path.append("..")
from multitaper_spectrogram_python import multitaper_spectrogram


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
    librosa.display.specshow(nanpow2db(root.spect), x_axis='time', y_axis='linear',
                             x_coords=root.stimes, y_coords=root.sfreqs, shading='auto', cmap="jet")
    plt.colorbar(label='Power (dB)')
    plt.xlabel("Time (HH:MM:SS)")
    plt.ylabel("Frequency (Hz)")
    plt.clim(clim)  # actually change colorbar scale
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
    multiprocess = False
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
    spath = resource_path("multitpaer_image.png")
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


if __name__ == '__main__':
    main()


# TODO: Button to export MT spect data, loading bar for ingesting edfs, loading icon/bar for spectrogram calc
