/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_multitaper_spectrogram_coder_mex.c
 *
 * Code generation for function '_coder_multitaper_spectrogram_coder_mex'
 *
 */

/* Include files */
#include "_coder_multitaper_spectrogram_coder_mex.h"
#include "_coder_multitaper_spectrogram_coder_api.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_initialize.h"
#include "multitaper_spectrogram_coder_terminate.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  static jmp_buf emlrtJBEnviron;
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexAtExit(&multitaper_spectrogram_coder_atexit);
  /* Initialize the memory manager. */
  omp_init_lock(&emlrtLockGlobal);
  omp_init_nest_lock(&multitaper_spectrogram_coder_nestLockGlobal);
  multitaper_spectrogram_coder_initialize();
  st.tls = emlrtRootTLSGlobal;
  emlrtSetJmpBuf(&st, &emlrtJBEnviron);
  if (setjmp(emlrtJBEnviron) == 0) {
    multitaper_spectrogram_coder_mexFunction(nlhs, plhs, nrhs, prhs);
    multitaper_spectrogram_coder_terminate();
    omp_destroy_lock(&emlrtLockGlobal);
    omp_destroy_nest_lock(&multitaper_spectrogram_coder_nestLockGlobal);
  } else {
    omp_destroy_lock(&emlrtLockGlobal);
    omp_destroy_nest_lock(&multitaper_spectrogram_coder_nestLockGlobal);
    emlrtReportParallelRunTimeError(&st);
  }
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal,
                           &emlrtLockerFunction, omp_get_num_procs(), NULL,
                           "UTF-8", true);
  return emlrtRootTLSGlobal;
}

void multitaper_spectrogram_coder_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                              int32_T nrhs,
                                              const mxArray *prhs[9])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[3];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 9) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 9, 4,
                        28, "multitaper_spectrogram_coder");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 28,
                        "multitaper_spectrogram_coder");
  }
  /* Call the function. */
  c_multitaper_spectrogram_coder_(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_multitaper_spectrogram_coder_mex.c) */
