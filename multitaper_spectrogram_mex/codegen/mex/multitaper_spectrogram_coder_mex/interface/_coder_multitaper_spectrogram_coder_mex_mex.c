/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_multitaper_spectrogram_coder_mex_mex.c
 *
 * Code generation for function '_coder_multitaper_spectrogram_coder_mex_mex'
 *
 */

/* Include files */
#include "_coder_multitaper_spectrogram_coder_mex_mex.h"
#include "_coder_multitaper_spectrogram_coder_mex_api.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex_initialize.h"
#include "multitaper_spectrogram_coder_mex_mexutil.h"
#include "multitaper_spectrogram_coder_mex_terminate.h"

/* Variable Definitions */
static jmp_buf emlrtJBEnviron;

/* Function Declarations */
MEXFUNCTION_LINKAGE void c_multitaper_spectrogram_coder_(int32_T nlhs, mxArray
  *plhs[3], int32_T nrhs, const mxArray *prhs[10]);

/* Function Definitions */
void c_multitaper_spectrogram_coder_(int32_T nlhs, mxArray *plhs[3], int32_T
  nrhs, const mxArray *prhs[10])
{
  const mxArray *outputs[3];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 10) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 10, 4,
                        32, "multitaper_spectrogram_coder_mex");
  }

  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 32,
                        "multitaper_spectrogram_coder_mex");
  }

  /* Call the function. */
  multitaper_spectrogram_coder_mex_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexAtExit(multitaper_spectrogram_coder_mex_atexit);
  emlrtLoadMATLABLibrary("sys/os/glnxa64/libiomp5.so");

  /* Initialize the memory manager. */
  omp_init_lock(&emlrtLockGlobal);
  omp_init_nest_lock(&emlrtNestLockGlobal);

  /* Module initialization. */
  multitaper_spectrogram_coder_mex_initialize();
  st.tls = emlrtRootTLSGlobal;
  emlrtSetJmpBuf(&st, &emlrtJBEnviron);
  if (setjmp(emlrtJBEnviron) == 0) {
    /* Dispatch the entry-point. */
    c_multitaper_spectrogram_coder_(nlhs, plhs, nrhs, prhs);

    /* Module termination. */
    multitaper_spectrogram_coder_mex_terminate();
    omp_destroy_lock(&emlrtLockGlobal);
    omp_destroy_nest_lock(&emlrtNestLockGlobal);
  } else {
    omp_destroy_lock(&emlrtLockGlobal);
    omp_destroy_nest_lock(&emlrtNestLockGlobal);
    emlrtReportParallelRunTimeError(&st);
  }
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal,
                     emlrtLockerFunction, omp_get_num_procs());
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_multitaper_spectrogram_coder_mex_mex.c) */
