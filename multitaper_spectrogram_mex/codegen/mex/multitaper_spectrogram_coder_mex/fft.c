/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fft.c
 *
 * Code generation for function 'fft'
 *
 */

/* Include files */
#include "fft.h"
#include "FFTWApi.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex_emxutil.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo me_emlrtRSI = { 8,  /* lineNo */
  "fft",                               /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/fft.m"/* pathName */
};

static emlrtRSInfo ne_emlrtRSI = { 102,/* lineNo */
  "fft",                               /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/fft.m"/* pathName */
};

static emlrtRSInfo oe_emlrtRSI = { 86, /* lineNo */
  "fft",                               /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/fft.m"/* pathName */
};

static emlrtRTEInfo o_emlrtRTEI = { 52,/* lineNo */
  27,                                  /* colNo */
  "fft",                               /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/fft.m"/* pName */
};

static emlrtRTEInfo p_emlrtRTEI = { 48,/* lineNo */
  27,                                  /* colNo */
  "fft",                               /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/fft.m"/* pName */
};

static emlrtRTEInfo q_emlrtRTEI = { 35,/* lineNo */
  27,                                  /* colNo */
  "fft",                               /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/fft.m"/* pName */
};

static emlrtRTEInfo qc_emlrtRTEI = { 86,/* lineNo */
  9,                                   /* colNo */
  "fft",                               /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/fft.m"/* pName */
};

/* Function Definitions */
void fft(const emlrtStack *sp, const emxArray_real32_T *x, real_T varargin_1,
         emxArray_creal32_T *y)
{
  int32_T sz_idx_0_tmp;
  int32_T i;
  int32_T loop_ub;
  int32_T i1;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &me_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  if (x->size[0] == 1) {
    emlrtErrorWithMessageIdR2018a(&st, &q_emlrtRTEI,
      "Coder:toolbox:autoDimIncompatibility",
      "Coder:toolbox:autoDimIncompatibility", 0);
  }

  if ((!(varargin_1 >= 0.0)) || (!(varargin_1 == muDoubleScalarFloor(varargin_1))))
  {
    emlrtErrorWithMessageIdR2018a(&st, &p_emlrtRTEI,
      "MATLAB:fftfcn:lengthNotNonNegInt", "MATLAB:fftfcn:lengthNotNonNegInt", 0);
  }

  if (!(varargin_1 <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&st, &o_emlrtRTEI, "MATLAB:pmaxsize",
      "MATLAB:pmaxsize", 0);
  }

  if ((x->size[0] == 0) || (x->size[1] == 0) || ((int32_T)varargin_1 == 0)) {
    b_st.site = &oe_emlrtRSI;
    sz_idx_0_tmp = (int32_T)varargin_1;
    if (sz_idx_0_tmp < 0) {
      emlrtNonNegativeCheckR2012b(sz_idx_0_tmp, &d_emlrtDCI, &b_st);
    }

    if (x->size[1] < 0) {
      emlrtNonNegativeCheckR2012b(x->size[1], &d_emlrtDCI, &b_st);
    }

    i = y->size[0] * y->size[1];
    y->size[0] = sz_idx_0_tmp;
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal32_T(&b_st, y, i, &qc_emlrtRTEI);
    if (sz_idx_0_tmp > x->size[0]) {
      sz_idx_0_tmp = x->size[1];
      for (i = 0; i < sz_idx_0_tmp; i++) {
        loop_ub = y->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          y->data[i1 + y->size[0] * i].re = 0.0F;
          y->data[i1 + y->size[0] * i].im = 0.0F;
        }
      }
    }
  } else {
    b_st.site = &ne_emlrtRSI;
    FFTWApi_fft1d(&b_st, x, (int32_T)varargin_1, x->size[1], y);
  }
}

/* End of code generation (fft.c) */
