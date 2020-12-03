/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * FFTWApi.c
 *
 * Code generation for function 'FFTWApi'
 *
 */

/* Include files */
#include "FFTWApi.h"
#include "emlrt.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo qe_emlrtRSI = { 28, /* lineNo */
  "FFTWApi/fft1d",                     /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+fftw/FFTWApi.m"/* pathName */
};

static emlrtRSInfo re_emlrtRSI = { 31, /* lineNo */
  "FFTWApi/fft1d",                     /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+fftw/FFTWApi.m"/* pathName */
};

static emlrtRSInfo ue_emlrtRSI = { 21, /* lineNo */
  "MATLABFFTWCallback/fft1d",          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+fftw/MATLABFFTWCallback.m"/* pathName */
};

static emlrtRTEInfo rc_emlrtRTEI = { 21,/* lineNo */
  17,                                  /* colNo */
  "MATLABFFTWCallback",                /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+fftw/MATLABFFTWCallback.m"/* pName */
};

/* Function Definitions */
void FFTWApi_fft1d(const emlrtStack *sp, const emxArray_real32_T *data, int32_T
                   fftlen, int32_T nfft, emxArray_creal32_T *y)
{
  int32_T i;
  int32_T loop_ub;
  int32_T b_loop_ub;
  int32_T i1;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &qe_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  st.site = &re_emlrtRSI;
  if (emlrtIsInParallelRegion(&st)) {
    emlrtFFTWSetNumThreads(1);
  } else {
    emlrtFFTWSetNumThreads(12);
  }

  b_st.site = &ue_emlrtRSI;
  if (fftlen < 0) {
    emlrtNonNegativeCheckR2012b(fftlen, &d_emlrtDCI, &b_st);
  }

  if (data->size[1] < 0) {
    emlrtNonNegativeCheckR2012b(data->size[1], &d_emlrtDCI, &b_st);
  }

  i = y->size[0] * y->size[1];
  y->size[0] = fftlen;
  y->size[1] = data->size[1];
  emxEnsureCapacity_creal32_T(&b_st, y, i, &rc_emlrtRTEI);
  if (fftlen > data->size[0]) {
    loop_ub = data->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_loop_ub = y->size[0];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        y->data[i1 + y->size[0] * i].re = 0.0F;
        y->data[i1 + y->size[0] * i].im = 0.0F;
      }
    }
  }

  emlrtFFTWF_1D_R2C(&data->data[0], (real32_T *)&y->data[0], 1, fftlen,
                    data->size[0], nfft, -1);
}

/* End of code generation (FFTWApi.c) */
