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
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo ve_emlrtRSI = { 64, /* lineNo */
  "fltpower",                          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/ops/power.m"/* pathName */
};

static emlrtRTEInfo sc_emlrtRTEI = { 64,/* lineNo */
  5,                                   /* colNo */
  "power",                             /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/ops/power.m"/* pName */
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
  st.site = &x_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  b_st.site = &ve_emlrtRSI;
  unnamed_idx_0 = (uint32_T)a->size[0];
  unnamed_idx_1 = (uint32_T)a->size[1];
  nx = y->size[0] * y->size[1];
  y->size[0] = (int32_T)unnamed_idx_0;
  y->size[1] = (int32_T)unnamed_idx_1;
  emxEnsureCapacity_real32_T(&b_st, y, nx, &sc_emlrtRTEI);
  c_st.site = &oc_emlrtRSI;
  nx = (int32_T)unnamed_idx_0 * (int32_T)unnamed_idx_1;
  d_st.site = &pc_emlrtRSI;
  if ((1 <= nx) && (nx > 2147483646)) {
    e_st.site = &q_emlrtRSI;
    check_forloop_overflow_error(&e_st);
  }

  for (k = 0; k < nx; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

/* End of code generation (power.c) */
