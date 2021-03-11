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

    if desired_type == "float":
        try:
            out = float(entry)
        except ValueError:
            print(entry_name + " could not be converted to a float")
            out = 0

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


def go_multitaper_spectrogram(root):
    """

    :param root:
    :return:
    """
    # Get data and fs from extracted edf data #
    channel = root.channel.get()
    index_channel = [i for i in range(len(root.signal_headers)) if root.signal_headers[i]['label'] == channel][0]
    fs = root.signal_headers[index_channel]['sample_rate']
    data = root.signals[index_channel]

    # Get all multitaper params and convert to appropriate data types #
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


def channel_changed(event, root):
    """

    :param event:
    :param root:
    :return:
    """

    # Read edf file
    loading_label = Label(root, text="Loading in EDF channel...")
    loading_label.place(relx=0.25, rely=0.2, anchor=CENTER)
    root.update()
    root.signals, root.signal_headers, root.header = pyedflib.highlevel.read_edf(root.filename)
    loading_label.configure(text='Channel Loaded')

    # Populate other parameters #
    # TODO - create button for info on all parameters
    # Freq range param
    Label(root, text="Frequency Range of Analysis:").place(relx=0.08, rely=0.275, anchor=CENTER)
    root.freq_from_entered = ttk.Entry(root, width=10)
    root.freq_from_entered.place(relx=0.2, rely=0.275, anchor=CENTER)
    root.freq_from_entered.delete(0, END)
    root.freq_from_entered.insert(0, 0)
    Label(root, text=" - ").place(relx=0.25, rely=0.275, anchor=CENTER)
    root.freq_to_entered = ttk.Entry(root, width=10)
    root.freq_to_entered.place(relx=0.3, rely=0.275, anchor=CENTER)
    root.freq_to_entered.delete(0, END)
    root.freq_to_entered.insert(0, 50)
    Label(root, text=" (Hz)").place(relx=0.35, rely=0.275, anchor=CENTER)

    # Taper Parameters
    Label(root, text="Half-Time Bandwidth (HW):").place(relx=0.0725, rely=0.31, anchor=CENTER)
    root.hw_entered = ttk.Entry(root, width=10)
    root.hw_entered.delete(0, END)
    root.hw_entered.insert(0, 5)
    root.hw_entered.place(relx=0.2, rely=0.31, anchor=CENTER)

    Label(root, text="Number of Tapers (L):").place(relx=0.0525, rely=0.35, anchor=CENTER)
    root.num_tapers_entered = ttk.Entry(root, width=10)
    root.num_tapers_entered.place(relx=0.15, rely=0.35, anchor=CENTER)
    root.num_tapers_entered.delete(0, END)
    root.num_tapers_entered.insert(0, 9)

    # Window Parameters
    y = 0.385
    Label(root, text="Window Size (N):").place(relx=0.04, rely=y, anchor=CENTER)
    root.w_size_entered = ttk.Entry(root, width=10)
    root.w_size_entered.place(relx=0.125, rely=y, anchor=CENTER)
    root.w_size_entered.delete(0, END)
    root.w_size_entered.insert(0, 10)
    Label(root, text="(s)").place(relx=0.15, rely=y, anchor=CENTER)

    y = 0.42
    Label(root, text="Window Step Size:").place(relx=0.0525, rely=y, anchor=CENTER)
    root.w_step_entered = ttk.Entry(root, width=10)
    root.w_step_entered.place(relx=0.14, rely=y, anchor=CENTER)
    root.w_step_entered.delete(0, END)
    root.w_step_entered.insert(0, 5)
    Label(root, text="(s)").place(relx=0.185, rely=y, anchor=CENTER)

    # Detrend method param
    y = 0.45
    Label(root, text="Window Detrend Method: ").place(relx=0.075, rely=y, anchor=CENTER)
    root.detrend_method_entered = ttk.Entry(root, width=10)
    root.detrend_method_entered.place(relx=0.175, rely=y, anchor=CENTER)
    root.detrend_method_entered.delete(0, END)
    root.detrend_method_entered.insert(0, 'linear')
    Label(root, text="Choices: 'linear' - default, 'constant', 'off'").place(relx=0.3, rely=y, anchor=CENTER)

    # Minimum NFFT Param
    y = 0.50
    Label(root, text="Minimum NFFT:").place(relx=0.075, rely=y, anchor=CENTER)
    root.min_nfft_entered = ttk.Entry(root, width=10)
    root.min_nfft_entered.place(relx=0.15, rely=y, anchor=CENTER)
    root.min_nfft_entered.delete(0, END)
    root.min_nfft_entered.insert(0, 0)

    # Taper Weighting Param
    y = 0.55
    Label(root, text="Taper Weighting Method:").place(relx=0.075, rely=y, anchor=CENTER)
    root.weighting_entered = ttk.Entry(root, width=10)
    root.weighting_entered.place(relx=0.175, rely=y, anchor=CENTER)
    root.weighting_entered.delete(0, END)
    root.weighting_entered.insert(0, 'unity')
    Label(root, text="Choices: 'unity' - default/recommended, 'eigen', 'adapt'").place(relx=0.35, rely=y, anchor=CENTER)

    # Multitaper Calc Button
    multitaper_go = Button(root, text="Calculate Multitaper Spectrogram", width=25, height=2, relief=GROOVE,
                           bg='#0064BB', fg="white", highlightbackground='#0064BB',
                           command=lambda: go_multitaper_spectrogram(root))
    multitaper_go.place(relx=0.5, rely=0.8, anchor=CENTER)


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
    label_channel.place(relx=0.05, rely=0.2, anchor=CENTER)

    # Create select channel dropdown
    selected_channel = StringVar()
    root.channel = ttk.Combobox(root, textvariable=selected_channel)
    root.channel['values'] = root.signal_labels
    root.channel['state'] = 'readonly'
    root.channel.place(relx=0.15, rely=0.2, anchor=CENTER)
    root.channel.bind('<<ComboboxSelected>>', lambda event: channel_changed(event, root))

    return root.filename


def main():
    # Create window
    root = Tk()
    root.geometry("1250x600")
    root.title('Multitaper Spectrogram Analysis')
    center(root)

    # Put image in window
    #spath = resource_path("Splash.png")
    #simg = ImageTk.PhotoImage(Image.open(spath))
    #panel = Label(root, image=simg)
    #center(root)
    #panel.pack()

    # Put button on window
    btn = Button(root, text="Choose EDF File", width=25, height=2, relief=GROOVE, bg='#0064BB',
                 fg="white", command=lambda: start_edf(root), highlightbackground='#0064BB')
    btn.place(relx=0.1, rely=0.1, anchor=CENTER)

    choose_file = Label(root, text="Select EDF file to begin")
    choose_file.place(relx=0.1, rely=0.05, anchor=CENTER)

    root.mainloop()


if __name__ == '__main__':
    main()