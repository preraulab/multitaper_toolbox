/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder_mex_mexutil.c
 *
 * Code generation for function 'multitaper_spectrogram_coder_mex_mexutil'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder_mex_mexutil.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
emlrtCTX emlrtGetRootTLSGlobal(void)
{
  return emlrtRootTLSGlobal;
}

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, const emlrtConstCTX aTLS,
  void *aData)
{
  omp_set_lock(&emlrtLockGlobal);
  emlrtCallLockeeFunction(aLockee, aTLS, aData);
  omp_unset_lock(&emlrtLockGlobal);
}

/* End of code generation (multitaper_spectrogram_coder_mex_mexutil.c) */
