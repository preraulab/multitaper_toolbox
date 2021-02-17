/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mean.c
 *
 * Code generation for function 'mean'
 *
 */

/* Include files */
#include "mean.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo se_emlrtRSI = { 49, /* lineNo */
  "mean",                              /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/mean.m"/* pathName */
};

/* Function Definitions */
void mean(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  int32_T vstride;
  int32_T xj;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &se_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  b_st.site = &te_emlrtRSI;
  if (x->size[0] == 0) {
    y->size[0] = 0;
  } else {
    c_st.site = &ue_emlrtRSI;
    vstride = x->size[0];
    xj = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&c_st, y, xj, &jd_emlrtRTEI);
    d_st.site = &ve_emlrtRSI;
    if (x->size[0] > 2147483646) {
      e_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    for (xj = 0; xj < vstride; xj++) {
      y->data[xj] = x->data[xj];
    }

    d_st.site = &we_emlrtRSI;
    if (x->size[0] > 2147483646) {
      e_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    for (xj = 0; xj < vstride; xj++) {
      y->data[xj] += x->data[vstride + xj];
    }
  }

  vstride = y->size[0];
  for (xj = 0; xj < vstride; xj++) {
    y->data[xj] /= 2.0F;
  }
}

/* End of code generation (mean.c) */
