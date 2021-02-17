/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * all.cpp
 *
 * Code generation for function 'all'
 *
 */

/* Include files */
#include "all.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo mb_emlrtRSI = { 13, /* lineNo */
  "all",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/all.m"/* pathName */
};

/* Function Definitions */
boolean_T all(const emlrtStack *sp, const emxArray_boolean_T *x)
{
  boolean_T y;
  int32_T ix;
  boolean_T exitg1;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &mb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  y = true;
  b_st.site = &nb_emlrtRSI;
  if ((1 <= x->size[1]) && (x->size[1] > 2147483646)) {
    c_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  ix = 1;
  exitg1 = false;
  while ((!exitg1) && (ix <= x->size[1])) {
    if (!x->data[ix - 1]) {
      y = false;
      exitg1 = true;
    } else {
      ix++;
    }
  }

  return y;
}

/* End of code generation (all.cpp) */
