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
static emlrtRSInfo
    od_emlrtRSI =
        {
            112,    /* lineNo */
            "mean", /* fcnName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
            "mean.m" /* pathName */
};

/* Function Definitions */
void mean(const emlrtStack *sp, const emxArray_real32_T *x,
          emxArray_real32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  int32_T vstride_tmp;
  int32_T xj;
  const real32_T *x_data;
  real32_T *y_data;
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
  y_data = y->data;
  x_data = x->data;
  st.site = &od_emlrtRSI;
  b_st.site = &pd_emlrtRSI;
  if (x->size[0] == 0) {
    y->size[0] = 0;
  } else {
    c_st.site = &qd_emlrtRSI;
    vstride_tmp = x->size[0];
    xj = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&c_st, y, xj, &fd_emlrtRTEI);
    y_data = y->data;
    d_st.site = &rd_emlrtRSI;
    if (x->size[0] > 2147483646) {
      e_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }
    for (xj = 0; xj < vstride_tmp; xj++) {
      y_data[xj] = x_data[xj];
    }
    d_st.site = &sd_emlrtRSI;
    for (xj = 0; xj < vstride_tmp; xj++) {
      y_data[xj] += x_data[vstride_tmp + xj];
    }
  }
  vstride_tmp = y->size[0];
  for (xj = 0; xj < vstride_tmp; xj++) {
    y_data[xj] /= 2.0F;
  }
}

/* End of code generation (mean.c) */
