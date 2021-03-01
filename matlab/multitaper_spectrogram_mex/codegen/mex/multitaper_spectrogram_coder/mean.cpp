/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mean.cpp
 *
 * Code generation for function 'mean'
 *
 */

/* Include files */
#include "mean.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo gf_emlrtRSI = { 49, /* lineNo */
  "mean",                              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/mean.m"/* pathName */
};

/* Function Definitions */
void mean(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T *y)
{
  int32_T vstride;
  int32_T i;
  int32_T xj;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &gf_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  b_st.site = &hf_emlrtRSI;
  if (x->size[0] == 0) {
    y->size[0] = 0;
  } else {
    c_st.site = &if_emlrtRSI;
    vstride = x->size[0];
    i = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&c_st, y, i, &de_emlrtRTEI);
    d_st.site = &jf_emlrtRSI;
    if (x->size[0] > 2147483646) {
      e_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    for (xj = 0; xj < vstride; xj++) {
      y->data[xj] = x->data[xj];
    }

    d_st.site = &kf_emlrtRSI;
    if (x->size[0] > 2147483646) {
      e_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    for (xj = 0; xj < vstride; xj++) {
      y->data[xj] += x->data[vstride + xj];
    }
  }

  xj = y->size[0];
  for (i = 0; i < xj; i++) {
    y->data[i] /= 2.0F;
  }
}

/* End of code generation (mean.cpp) */
