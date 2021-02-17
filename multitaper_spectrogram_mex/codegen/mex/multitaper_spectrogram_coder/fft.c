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
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo je_emlrtRSI = { 63, /* lineNo */
  "fft",                               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/fft.m"/* pathName */
};

static emlrtRSInfo ke_emlrtRSI = { 31, /* lineNo */
  "fft",                               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+fft/fft.m"/* pathName */
};

static emlrtRSInfo le_emlrtRSI = { 58, /* lineNo */
  "executeCallback",                   /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+fft/fft.m"/* pathName */
};

static emlrtRSInfo me_emlrtRSI = { 44, /* lineNo */
  "Custom1DFFTCallback/fft",           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/shared/coder/coder/lib/+coder/+internal/Custom1DFFTCallback.m"/* pathName */
};

static emlrtRSInfo ne_emlrtRSI = { 54, /* lineNo */
  "Custom1DFFTCallback/fftLoop",       /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/shared/coder/coder/lib/+coder/+internal/Custom1DFFTCallback.m"/* pathName */
};

static emlrtRTEInfo m_emlrtRTEI = { 52,/* lineNo */
  35,                                  /* colNo */
  "fft",                               /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/fft.m"/* pName */
};

static emlrtRTEInfo n_emlrtRTEI = { 48,/* lineNo */
  35,                                  /* colNo */
  "fft",                               /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/fft.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 37,/* lineNo */
  31,                                  /* colNo */
  "fft",                               /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/fft.m"/* pName */
};

static emlrtRTEInfo gd_emlrtRTEI = { 63,/* lineNo */
  5,                                   /* colNo */
  "fft",                               /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/fft.m"/* pName */
};

static emlrtRTEInfo hd_emlrtRTEI = { 26,/* lineNo */
  32,                                  /* colNo */
  "MATLABFFTWCallback",                /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+fft/MATLABFFTWCallback.m"/* pName */
};

/* Function Definitions */
void fft(const emlrtStack *sp, const emxArray_real32_T *x, real_T varargin_1,
         emxArray_creal32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  int32_T i;
  int32_T loop_ub;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  if (x->size[0] == 1) {
    emlrtErrorWithMessageIdR2018a(sp, &o_emlrtRTEI,
      "Coder:toolbox:autoDimIncompatibility",
      "Coder:toolbox:autoDimIncompatibility", 0);
  }

  if (!(varargin_1 == muDoubleScalarFloor(varargin_1))) {
    emlrtErrorWithMessageIdR2018a(sp, &n_emlrtRTEI,
      "MATLAB:fftfcn:lengthNotNonNegInt", "MATLAB:fftfcn:lengthNotNonNegInt", 0);
  }

  if (!(varargin_1 <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(sp, &m_emlrtRTEI, "MATLAB:pmaxsize",
      "MATLAB:pmaxsize", 0);
  }

  st.site = &je_emlrtRSI;
  if ((x->size[0] == 0) || (x->size[1] == 0) || (0 == (int32_T)varargin_1)) {
    i = y->size[0] * y->size[1];
    y->size[0] = (int32_T)varargin_1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal32_T(&st, y, i, &gd_emlrtRTEI);
    loop_ub = (int32_T)varargin_1 * x->size[1];
    for (i = 0; i < loop_ub; i++) {
      y->data[i].re = 0.0F;
      y->data[i].im = 0.0F;
    }
  } else {
    b_st.site = &ke_emlrtRSI;
    c_st.site = &le_emlrtRSI;
    d_st.site = &me_emlrtRSI;
    e_st.site = &ne_emlrtRSI;
    if (emlrtIsInParallelRegion(&e_st)) {
      emlrtFFTWSetNumThreads(1);
    } else {
      emlrtFFTWSetNumThreads(28);
    }

    i = y->size[0] * y->size[1];
    y->size[0] = (int32_T)varargin_1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal32_T(&e_st, y, i, &hd_emlrtRTEI);
    emlrtFFTWF_1D_R2C(&x->data[0], (real32_T *)&y->data[0], 1, (int32_T)
                      varargin_1, x->size[0], x->size[1], -1);
  }
}

/* End of code generation (fft.c) */
