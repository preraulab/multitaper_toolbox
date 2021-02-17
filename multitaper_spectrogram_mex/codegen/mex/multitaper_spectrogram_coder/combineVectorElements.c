/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * combineVectorElements.c
 *
 * Code generation for function 'combineVectorElements'
 *
 */

/* Include files */
#include "combineVectorElements.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo kb_emlrtRSI = { 133,/* lineNo */
  "combineVectorElements",             /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

static emlrtRSInfo lb_emlrtRSI = { 194,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

/* Function Definitions */
int32_T combineVectorElements(const emlrtStack *sp, const emxArray_boolean_T *x)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  int32_T k;
  int32_T vlen;
  int32_T y;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  vlen = x->size[1];
  if (x->size[1] == 0) {
    y = 0;
  } else {
    st.site = &kb_emlrtRSI;
    y = x->data[0];
    b_st.site = &lb_emlrtRSI;
    if ((2 <= x->size[1]) && (x->size[1] > 2147483646)) {
      c_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (k = 2; k <= vlen; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

/* End of code generation (combineVectorElements.c) */
