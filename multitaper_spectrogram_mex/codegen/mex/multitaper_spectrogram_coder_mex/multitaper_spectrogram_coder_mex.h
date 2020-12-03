/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder_mex.h
 *
 * Code generation for function 'multitaper_spectrogram_coder_mex'
 *
 */

#pragma once

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "omp.h"
#include "multitaper_spectrogram_coder_mex_types.h"

/* Function Declarations */
void multitaper_spectrogram_coder_mex(const emlrtStack *sp, const
  emxArray_real32_T *data, real_T Fs, real_T frequency_range[2], const
  emxArray_real_T *DPSS_tapers, real_T time_bandwidth, real_T num_tapers, real_T
  winsize_samples, real_T winstep_samples, real_T min_NFFT, real_T detrend_opt,
  emxArray_real32_T *mt_spectrogram, emxArray_real_T *stimes, emxArray_real_T
  *sfreqs);

/* End of code generation (multitaper_spectrogram_coder_mex.h) */
