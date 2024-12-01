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
static emlrtRSInfo fd_emlrtRSI = {
    63,    /* lineNo */
    "fft", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/fft.m" /* pathName
                                                                            */
};

static emlrtRSInfo gd_emlrtRSI = {
    42,    /* lineNo */
    "fft", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+fft/"
    "fft.m" /* pathName */
};

static emlrtRSInfo hd_emlrtRSI = {
    69,                /* lineNo */
    "executeCallback", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+fft/"
    "fft.m" /* pathName */
};

static emlrtRSInfo id_emlrtRSI = {
    44,                        /* lineNo */
    "Custom1DFFTCallback/fft", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/shared/coder/coder/lib/+coder/"
    "+internal/Custom1DFFTCallback.m" /* pathName */
};

static emlrtRSInfo jd_emlrtRSI = {
    54,                            /* lineNo */
    "Custom1DFFTCallback/fftLoop", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/shared/coder/coder/lib/+coder/"
    "+internal/Custom1DFFTCallback.m" /* pathName */
};

static emlrtRTEInfo i_emlrtRTEI = {
    52,    /* lineNo */
    35,    /* colNo */
    "fft", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/fft.m" /* pName
                                                                            */
};

static emlrtRTEInfo j_emlrtRTEI = {
    48,    /* lineNo */
    35,    /* colNo */
    "fft", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/fft.m" /* pName
                                                                            */
};

static emlrtRTEInfo k_emlrtRTEI = {
    37,    /* lineNo */
    31,    /* colNo */
    "fft", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/fft.m" /* pName
                                                                            */
};

static emlrtRTEInfo dd_emlrtRTEI = {
    63,    /* lineNo */
    5,     /* colNo */
    "fft", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/fft.m" /* pName
                                                                            */
};

static emlrtRTEInfo ed_emlrtRTEI = {
    32,                   /* lineNo */
    36,                   /* colNo */
    "MATLABFFTWCallback", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+fft/"
    "MATLABFFTWCallback.m" /* pName */
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
  creal32_T *y_data;
  int32_T i;
  const real32_T *x_data;
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
  x_data = x->data;
  if (x->size[0] == 1) {
    emlrtErrorWithMessageIdR2018a(sp, &k_emlrtRTEI,
                                  "Coder:toolbox:autoDimIncompatibility",
                                  "Coder:toolbox:autoDimIncompatibility", 0);
  }
  if (!(varargin_1 == muDoubleScalarFloor(varargin_1))) {
    emlrtErrorWithMessageIdR2018a(sp, &j_emlrtRTEI,
                                  "MATLAB:fftfcn:lengthNotNonNegInt",
                                  "MATLAB:fftfcn:lengthNotNonNegInt", 0);
  }
  if (!(varargin_1 <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(sp, &i_emlrtRTEI, "MATLAB:pmaxsize",
                                  "MATLAB:pmaxsize", 0);
  }
  st.site = &fd_emlrtRSI;
  if ((x->size[0] == 0) || (x->size[1] == 0) || ((int32_T)varargin_1 == 0)) {
    int32_T loop_ub_tmp;
    i = y->size[0] * y->size[1];
    y->size[0] = (int32_T)varargin_1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal32_T(&st, y, i, &dd_emlrtRTEI);
    y_data = y->data;
    loop_ub_tmp = (int32_T)varargin_1 * x->size[1];
    for (i = 0; i < loop_ub_tmp; i++) {
      y_data[i].re = 0.0F;
      y_data[i].im = 0.0F;
    }
  } else {
    b_st.site = &gd_emlrtRSI;
    c_st.site = &hd_emlrtRSI;
    d_st.site = &id_emlrtRSI;
    e_st.site = &jd_emlrtRSI;
    if (emlrtIsInParallelRegion(&e_st)) {
      emlrtFFTWSetNumThreads(1);
    } else {
      emlrtFFTWSetNumThreads(16);
    }
    i = y->size[0] * y->size[1];
    y->size[0] = (int32_T)varargin_1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal32_T(&e_st, y, i, &ed_emlrtRTEI);
    y_data = y->data;
    emlrtFFTWF_1D_R2C((real32_T *)&x_data[0], (real32_T *)&y_data[0], 1,
                      (int32_T)varargin_1, x->size[0], x->size[1], -1);
  }
}

/* End of code generation (fft.c) */
