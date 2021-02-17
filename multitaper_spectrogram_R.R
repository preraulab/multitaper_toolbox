## Package installation and library initializing ##
packages_req <- c("waveslim", "pracma", "fields")
not_installed <- packages_req[!(packages_req %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)                               # Install not installed packages

library(waveslim)
library(pracma)
library(fields)


# Multitaper Spectrogram #
multitaper_spectrogram_R <- function(data, fs, frequency_range=NULL, time_bandwidth=5, num_tapers=NULL, window_params=c(5,1),
                                     min_nfft=0, detrend_opt='linear', plot_on=TRUE, verbose=TRUE){
  # Compute multitaper spectrogram of timeseries data
  # 
  # Results tend to agree with Prerau Lab python implementation of multitaper spectrogram with precision on the order of at most 
  # 10^-7 with SD of at most 10^-5
  #
  # params:
  #         data (numeric vector): time series data -- required
  #         fs (numeric): sampling frequency in Hz  -- required
  #         frequency_range (numeric vector): c(<min frequency>, <max frequency>) (default: NULL, adjusted to 
  #                                           c(0, nyquist) later)
  #         time_bandwidth (numeric): time-half bandwidth product (window duration*half bandwidth of main lobe)
  #                                   (default: 5 Hz*s)
  #         num_tapers (numeric): number of DPSS tapers to use (default: NULL [will be computed
  #                                                               as floor(2*time_bandwidth - 1)])
  #         window_params (numeric vector): c(window size (seconds), step size (seconds)) (default: [5 1])
  #         detrend_opt (char): detrend data window ('linear' (default), 'constant', 'off') 
  #         min_nfft (numeric): minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0)
  #         plot_on (logical): plot results (default: TRUE)
  #         verbose (logical): display spectrogram properties (default: TRUE)
  #
  # returns:
  #         mt_spectrogram (matrix): spectral power matrix
  #         stimes (numeric vector): timepoints (s) in mt_spectrogram
  #         sfreqs (numeric vector): frequency values (Hz) in mt_spectrogram
  
  # Process user input
  res <- process_input(data, fs, frequency_range, time_bandwidth, num_tapers, window_params, min_nfft, detrend_opt, 
                       plot_on, verbose)
  
  data <- res[[1]]
  fs <- res[[2]]
  frequency_range <- res[[3]]
  time_bandwidth <- res[[4]]
  num_tapers <- res[[5]]
  winsize_samples <- res[[6]]
  winstep_samples = res[[7]]
  window_start = res[[8]]
  num_windows <- res[[9]]
  nfft <- res[[10]]
  detrend_opt <- res[[11]]
  plot_on <- res[[12]]
  verbose <- res[[13]]
  
  # Set up spectrogram parameters
  res <- process_spectrogram_params(fs, nfft, frequency_range, window_start, winsize_samples)
  window_idxs <- res[[1]]
  stimes <- res[[2]]
  sfreqs <- res[[3]]
  freq_inds <- res[[4]]
  
  # Display spectrogram parameters if desired
  if(verbose){
    display_spectrogram_properties(fs, time_bandwidth, num_tapers, c(winsize_samples, winstep_samples), frequency_range, 
                                   detrend_opt)
  }
  
  # Split data into window segments
  data_segments <- t(sapply(window_idxs, split_data_helper, data=data))
  
  # COMPUTE THE MULTITAPER SPECTROGRAM
  #     STEP 1: Compute DPSS tapers based on desired spectral properties
  #     STEP 2: Multiply the data segment by the DPSS Tapers
  #     STEP 3: Compute the spectrum for each tapered segment
  #     STEP 4: Take the mean of the tapered spectra
  
  # Compute DPSS tapers (STEP 1)
  dpss_tapers <- dpss.taper(winsize_samples, num_tapers, time_bandwidth) * sqrt(fs)
  
  tic <- proc.time() # start timer for multitaper
  # Compute multitaper
  mt_spectrogram = apply(data_segments, 1, calc_mts_segment, dpss_tapers=dpss_tapers, nfft=nfft, 
                         freq_inds=freq_inds, detrend_opt=detrend_opt)
  
  # Compute mean fft magnitude (STEP 4.2)
  mt_spectrogram = Conj(t(mt_spectrogram)) / fs^2 / num_tapers
  
  # End timer and get elapsed time
  toc = proc.time()
  elapsed = toc-tic
  if(verbose){
    print(paste("Multitaper compute time: ", toString(round(elapsed[[3]], digits=5)), " seconds", sep=""))
  }
  
  
  if(all(as.vector(mt_spectrogram) == 0)){
    print("Spectrogram calculated as all zeros, no plot shown")
  }else if(plot_on){
    image.plot(x=stimes, y=sfreqs, nanpow2db(mt_spectrogram), xlab="Time (s)", 
               ylab='Frequency (Hz)')
  }
  
  return(list(mt_spectrogram, stimes, sfreqs))
}

split_data_helper <- function(indices, data){ # for sapply when splitting data into windows
  data_seg = data[indices]
  return(data_seg)
}



### Helper Functions ###

# Process user input #
process_input <- function(data, fs, frequency_range=NULL, time_bandwidth=5, num_tapers=NULL,
                          window_params=c(5,1), min_nfft=0, detrend_opt='linear', plot_on=TRUE,
                          verbose=TRUE){
  
  # Helper function to process multitaper_spectrogram arguments, mainly checking for validity
  #
  # Params:
  #        data (numeric vector): time series data -- required
  #        fs (numeric): sampling frequency in Hz -- required
  #        frequency range (numeric vector): c(<min frequency>, <max frequency>) (default: c(0 nyquist))
  #        time_bandwidth (numeric): time-half bandwidth product (window duration*half bandwidth of main lobe) (default: 5 Hz*s)
  #        num_tapers (numeric): number of DPSS tapers to use (default None [will be computed as floor(2*time_bandwidth - 1)])
  #        window_params (numeric vector): c(window size (seconds), step size (seconds)) default: c(5,1)
  #        detrend_opt (char): detrend data window ('linear' (default), 'constant', 'off')
  #        min_nfft: (numeric): minimum allowable NFFT size, adds zero padding for interpolation (default: 0)
  #        plot_on: (logical): plot results (default: TRUE)
  #        verbose (logical)L display spectrogram properties (default; TRUE)
  #
  # Returns:
  #         data (numeric vector) same as input
  #         fs (numeric): same as input
  #         frequency_range (numeric vector): same as input or calculated from fs if not given
  #         time_bandwidth (numeric): same as input or default if not given
  #         num_tapers (numeric): same as input or calculated from time bandwidth if not given
  #         winsize_samples (numeric): number of samples in a single time window
  #         winstep_samples (numeric): number of samples in a single window step
  #         window_start (numeric vector): matrix of timestamps representing the beginning time for each window
  #         num_windows (numeric): number of total windows
  #         nfft (numeric): length of signal to calculate fft on
  #         detrend_opt (char): same as input or default if not given
  #         plot_on (logical): same as input or default if not given
  #         verbose (logical): same as input or default if not given
  
  
  # Make sure data is 1D atomic vector
  if((is.atomic(data) == FALSE) | is.list(data)){
    stop("data must be a 1D atomic vector")
  }
  
  # Set frequency range if not provided
  if(is.null(frequency_range)){
    frequency_range <- c(0, fs/2)
  }
  
  # Set detrend method
  detrend_opt = tolower(detrend_opt)
  if(detrend_opt != 'linear'){
    if(detrend_opt == 'const'){
      detrend_opt <- 'constant'
    } else if(detrend_opt == 'none' || detrend_opt == 'false'){
      detrend_opt <- 'off'
    }else{
      stop(paste("'", toString(detrend_opt), "' is not a valid detrend_opt argument. The",
                 " choices are: 'constant', 'linear', or 'off'.", sep=""))
    } 
  }
  
  # Check if frequency range is valid
  if(frequency_range[2] > fs/2){
    frequency_range[2] <- fs/2
    warning(paste("Upper frequency range greater than Nyquist, setting range to [",
                  toString(frequency_range[1]), ",", toString(frequency_range[2]), "].",
                  sep=""))
  }
  
  # Set number of tapers if none provided
  optimal_num_tapers = floor(2*time_bandwidth) - 1
  if(is.null(num_tapers)){
    num_tapers <- optimal_num_tapers
  }
  
  # Warn if number of tapers is suboptimal
  if(num_tapers != optimal_num_tapers){
    warning(paste("Suboptimal number of tapers being used. Number of tapers is optimal at floor(2*TW) - 1 which is ",
                  toString(optimal_num_tapers), " in this case.", sep=""))
  }
  
  
  # Check if window size is valid, fix if not
  if((window_params[1]*fs) %% 1 != 0){
    winsize_samples <- round(window_params[1]*fs)
    warning(paste("Window size is not divisible by sampling frequency. Adjusting window",
                  " size to ", toString(winsize_samples/fs), " seconds.", sep=""))
  } else{
    winsize_samples <- window_params[1]*fs
  }
  
  # Check if window step size is valid, fix if not
  if((window_params[2]*fs) %% 1 != 0){
    winstep_samples <- round(window_params[2]*fs)
    warning(paste("Window step size is not divisible by sampling frequency. Adjusting window",
                  " step size to ", toString(winstep_samples/fs), " seconds.", sep=""))
  } else{
    winstep_samples <- window_params[2]*fs
  }
  
  # Get total data length
  len_data = length(data)
  
  # Check if length of data is smaller than window (bad)
  if(len_data < winsize_samples){
    stop(paste("Data length (", toString(len_data), ") is shorter than the window size (",
               toString(winsize_samples), "). Either increase data length or decrease",
               " window size.", sep=""))
  }
  
  # Find window start indices and num of windows
  window_start = seq(1, len_data-winsize_samples+1, by=winstep_samples)
  num_windows = length(window_start)
  
  # Get num points in FFT
  nfft = max(max(2^ceiling(log2(abs(winsize_samples))), winsize_samples), 2^ceiling(log2(abs(min_nfft))))
  
  return(list(data, fs, frequency_range, time_bandwidth, num_tapers, winsize_samples, winstep_samples, 
              window_start, num_windows, nfft, detrend_opt, plot_on, verbose))
}



# Process spectrogram inputs #
process_spectrogram_params <- function(fs, nfft, frequency_range, window_start, datawin_size){
  # Helper function to create frequency vector and window indices
  #
  # Params:
  #         fs (numeric): sampling frequency in Hz  -- required
  #         nfft (numeric): length of signal to calculate fft on -- required
  #         window_start (numeric vector): timestamps representing the beginning time for each window -- required
  #         datawin_size (numeric): seconds in one window -- required
  # 
  # Returns:
  #         window_idxs (matrix): indices of timestamps for each window (nxm where n=number of windows and m=datawin_size)
  #         stimes (numeric vector): times for the centers of the spectral bins (1xt)
  #         sfreqs (numeric vector): frequency bins for spectrogram (1xf)
  #         freq_inds (logical vector): indicates which frequencies are being analyzed in an array of frequencies from 0 to fs
  #                    with steps of fs/nfft
  
  
  # Create frequency vector
  df <- fs/nfft
  sfreqs <- seq(df/2, fs-(df/2), by=df)
  
  # Get frequencies for given frequency range
  freq_inds <- (sfreqs >= frequency_range[1]) & (sfreqs <= frequency_range[2])
  sfreqs <- sfreqs[freq_inds]
  
  # Compute times in middle of each spectrum
  window_middle_times <- window_start + round(datawin_size/2)
  stimes <- window_middle_times / fs
  
  # Get indices for each window
  window_idxs <- lapply(window_start, window_index_helper, datawin_size=datawin_size) # list of indices for n windows
  
  
  return(list(window_idxs, stimes, sfreqs, freq_inds))
  
}

window_index_helper <- function(start, datawin_size){
  res = seq(start, start+datawin_size-1, by=1)
  return(res)
}



# Display Spectrogram Properties #
display_spectrogram_properties <- function(fs, time_bandwidth, num_tapers, data_window_params, frequency_range, detrend_opt){
  # Prints spectrogram properties
  #
  # Params:
  #         fs (numeric): sampling frequency in Hz  -- required
  #         time_bandwidth (numeric): time-half bandwidth product (window duration*1/2*frequency_resolution) -- required
  #         num_tapers (numeric): number of DPSS tapers to use -- required
  #         data_window_params (numeric vector): c(window length(s), window step size(s) -- required
  #         frequency_range (numeric vector): c(<min frequency>, <max frequency>) -- required
  #         detrend_opt (char): detrend data window ('linear' (default), 'constant', 'off')
  #
  # Returns:
  #         This function does not return anythin
  
  data_window_params = data_window_params / fs
  
  # Print spectrogram properties
  print("Multitaper Spectrogram Properties: ")
  print(paste('     Spectral Resolution: ', toString(2 * time_bandwidth / data_window_params[1]), 'Hz', sep=""))
  print(paste('     Window Length: ', toString(data_window_params[1]), 's', sep=""))
  print(paste('     Window Step: ', toString(data_window_params[2]), 's', sep=""))
  print(paste('     Time Half-Bandwidth Product: ', toString(time_bandwidth), sep=""))
  print(paste('     Number of Tapers: ', toString(num_tapers), sep=""))
  print(paste('     Frequency Range: ', toString(frequency_range[1]), "-", toString(frequency_range[2]), 'Hz', sep=""))
  print(paste('     Detrend: ', detrend_opt, sep=""))
  
}


# Convert power to dB #
nanpow2db <- function(y){
  # Power to dB conversion, setting negatives and zeros to NaN
  #
  # params: 
  #         y: power --required
  #
  # returns:
  #         ydB: dB (with 0s and negativs set to NaN)
  
  if(length(y)==1){
    if(y==0){
      return(NaN)
    } else(ydB <- 10*log10(y))
  }else{
    y[y==0] <- NaN
    ydB <- 10*log10(y)
  }
  return(ydB)
}



# Calculate multitpaer spectrum of single segment #
calc_mts_segment <- function(data_segment, dpss_tapers, nfft, freq_inds, detrend_opt){
  # Calculate multitaper spectrum for a single segment of data
  #
  # params:
  #         data_segment (numeric vector): segment of the EEG data of length window size (s) * fs -- required
  #         dpss_tapers (numeric matrix): DPSS taper params to multiply signal by. Dims are (num_tapers, winsize_samples)
  #                                        -- required
  #         nfft (numeric): length of signal to calculate fft on -- required 
  #         freq_inds (logical vector): boolean array indicating frequencies to use in an array of frequenices
  #                                    from 0 to fs with steps of fs/nfft --required
  #         detrend_opt (char): detrend data window ('linear' (default), 'constant', 'off') --required
  #
  # returns:
  #         mt_spectrum (numeric matrix): spectral power for single window
  
  # If segment has all zeros, return vector of zeros
  if(all(data_segment==0)){
    ret <- rep(0, sum(freq_inds))
    return(ret)
  }
  
  # Optionally detrend data to remove low freq DC component
  if(detrend_opt != 'off'){
    data_segment <- detrend(data_segment, tt=detrend_opt)
  }
  
  # Multiply data by dpss tapers (STEP 2)
  tapered_data <- sweep(dpss_tapers, 1, data_segment, '*')
  
  # Manually add nfft zero-padding (R's fft function does not support)
  tapered_padded_data <- rbind(tapered_data, matrix(0, nrow=nfft-nrow(tapered_data), ncol=ncol(tapered_data)))
  
  # Compute the FFT (STEP 3)
  fft_data <- apply(tapered_padded_data, 2, fft)
  fft_range = fft_data[freq_inds,]
  
  # Take the FFT magnitude (STEP 4.1)
  magnitude = Im(fft_range)^2 + Re(fft_range)^2
  mt_spectrum = rowSums(magnitude)
  
  return(mt_spectrum)
}

