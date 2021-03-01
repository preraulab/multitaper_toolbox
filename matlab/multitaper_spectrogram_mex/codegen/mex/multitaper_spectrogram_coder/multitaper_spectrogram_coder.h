/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder.h
 *
 * Code generation for function 'multitaper_spectrogram_coder'
 *
 */

#pragma once

/* Include files */
#include "multitaper_spectrogram_coder_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include "omp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
emlrtCTX emlrtGetRootTLSGlobal(void);
void emlrtLockerFunction(EmlrtLockeeFunction aLockee, const emlrtConstCTX aTLS,
  void *aData);
void multitaper_spectrogram_coder(const emlrtStack *sp, const emxArray_real32_T *
  data, real_T Fs, const real_T frequency_range[2], const emxArray_real_T
  *DPSS_tapers, const emxArray_real_T *DPSS_eigen, real_T winstep_samples,
  real_T min_NFFT, real_T detrend_opt, real_T weighting, emxArray_real32_T
  *mt_spectrogram, emxArray_real_T *stimes, emxArray_real_T *sfreqs);

/* End of code generation (multitaper_spectrogram_coder.h) */
