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
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo bb_emlrtRSI = { 124,/* lineNo */
  "combineVectorElements",             /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

static emlrtRSInfo cb_emlrtRSI = { 184,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

/* Function Definitions */
int32_T combineVectorElements(const emlrtStack *sp, const emxArray_boolean_T *x)
{
  int32_T y;
  int32_T vlen;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
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
    st.site = &bb_emlrtRSI;
    y = x->data[0];
    b_st.site = &cb_emlrtRSI;
    if ((2 <= x->size[1]) && (x->size[1] > 2147483646)) {
      c_st.site = &q_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (k = 2; k <= vlen; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

/* End of code generation (combineVectorElements.c) */
