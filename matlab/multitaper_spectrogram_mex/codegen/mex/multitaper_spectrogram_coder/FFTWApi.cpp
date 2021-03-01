/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * FFTWApi.cpp
 *
 * Code generation for function 'FFTWApi'
 *
 */

/* Include files */
#include "FFTWApi.h"
#include "emlrt.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo ye_emlrtRSI = { 31, /* lineNo */
  "FFTWApi/fft1d",                     /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+fftw/FFTWApi.m"/* pathName */
};

static emlrtRTEInfo ae_emlrtRTEI = { 21,/* lineNo */
  17,                                  /* colNo */
  "MATLABFFTWCallback",                /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+fftw/MATLABFFTWCallback.m"/* pName */
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
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &ye_emlrtRSI;
  if (emlrtIsInParallelRegion(&st)) {
    emlrtFFTWSetNumThreads(1);
  } else {
    emlrtFFTWSetNumThreads(4);
  }

  i = y->size[0] * y->size[1];
  y->size[0] = fftlen;
  y->size[1] = data->size[1];
  emxEnsureCapacity_creal32_T(&st, y, i, &ae_emlrtRTEI);
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

/* End of code generation (FFTWApi.cpp) */
