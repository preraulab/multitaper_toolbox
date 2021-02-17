/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * power.c
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "power.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo re_emlrtRSI = { 79, /* lineNo */
  "fltpower",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/power.m"/* pathName */
};

/* Function Definitions */
void power(const emlrtStack *sp, const emxArray_real32_T *a, emxArray_real32_T
           *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  int32_T k;
  int32_T nx;
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
  b_st.site = &re_emlrtRSI;
  nx = y->size[0] * y->size[1];
  y->size[0] = a->size[0];
  y->size[1] = a->size[1];
  emxEnsureCapacity_real32_T(&b_st, y, nx, &id_emlrtRTEI);
  c_st.site = &xc_emlrtRSI;
  nx = a->size[0] * a->size[1];
  d_st.site = &yc_emlrtRSI;
  if ((1 <= nx) && (nx > 2147483646)) {
    e_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&e_st);
  }

  for (k = 0; k < nx; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

/* End of code generation (power.c) */
