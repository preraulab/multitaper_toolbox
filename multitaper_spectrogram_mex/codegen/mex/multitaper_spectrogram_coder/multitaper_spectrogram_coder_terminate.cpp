/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder_terminate.cpp
 *
 * Code generation for function 'multitaper_spectrogram_coder_terminate'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder_terminate.h"
#include "_coder_multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void multitaper_spectrogram_coder_atexit()
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void multitaper_spectrogram_coder_terminate()
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (multitaper_spectrogram_coder_terminate.cpp) */
