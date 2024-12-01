/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder_initialize.c
 *
 * Code generation for function 'multitaper_spectrogram_coder_initialize'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder_initialize.h"
#include "_coder_multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void multitaper_spectrogram_coder_once(void);

/* Function Definitions */
static void multitaper_spectrogram_coder_once(void)
{
  mex_InitInfAndNan();
}

void multitaper_spectrogram_coder_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    multitaper_spectrogram_coder_once();
  }
}

/* End of code generation (multitaper_spectrogram_coder_initialize.c) */
