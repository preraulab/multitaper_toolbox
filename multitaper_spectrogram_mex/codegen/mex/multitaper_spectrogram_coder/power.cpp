/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * power.cpp
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "power.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo cf_emlrtRSI = { 64, /* lineNo */
  "fltpower",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/power.m"/* pathName */
};

static emlrtRTEInfo be_emlrtRTEI = { 64,/* lineNo */
  5,                                   /* colNo */
  "power",                             /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/power.m"/* pName */
};

/* Function Definitions */
void power(const emlrtStack *sp, const emxArray_real32_T *a, emxArray_real32_T
           *y)
{
  uint32_T unnamed_idx_0;
  uint32_T unnamed_idx_1;
  int32_T nx;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &hb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  b_st.site = &cf_emlrtRSI;
  unnamed_idx_0 = static_cast<uint32_T>(a->size[0]);
  unnamed_idx_1 = static_cast<uint32_T>(a->size[1]);
  nx = y->size[0] * y->size[1];
  y->size[0] = static_cast<int32_T>(unnamed_idx_0);
  y->size[1] = static_cast<int32_T>(unnamed_idx_1);
  emxEnsureCapacity_real32_T(&b_st, y, nx, &be_emlrtRTEI);
  c_st.site = &wc_emlrtRSI;
  nx = static_cast<int32_T>(unnamed_idx_0) * static_cast<int32_T>(unnamed_idx_1);
  d_st.site = &xc_emlrtRSI;
  if ((1 <= nx) && (nx > 2147483646)) {
    e_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&e_st);
  }

  for (k = 0; k < nx; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

/* End of code generation (power.cpp) */
