/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find.c
 *
 * Code generation for function 'find'
 *
 */

/* Include files */
#include "find.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo hf_emlrtRSI = { 144,/* lineNo */
  "eml_find",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pathName */
};

static emlrtRSInfo if_emlrtRSI = { 382,/* lineNo */
  "find_first_indices",                /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pathName */
};

static emlrtRTEInfo r_emlrtRTEI = { 392,/* lineNo */
  1,                                   /* colNo */
  "find_first_indices",                /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pName */
};

static emlrtRTEInfo nd_emlrtRTEI = { 364,/* lineNo */
  24,                                  /* colNo */
  "find",                              /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pName */
};

static emlrtRTEInfo od_emlrtRTEI = { 144,/* lineNo */
  9,                                   /* colNo */
  "find",                              /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pName */
};

/* Function Definitions */
void eml_find(const emlrtStack *sp, const emxArray_boolean_T *x,
              emxArray_int32_T *i)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  int32_T idx;
  int32_T ii;
  int32_T nx;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  nx = x->size[1];
  st.site = &hf_emlrtRSI;
  idx = 0;
  ii = i->size[0] * i->size[1];
  i->size[0] = 1;
  i->size[1] = x->size[1];
  emxEnsureCapacity_int32_T(&st, i, ii, &nd_emlrtRTEI);
  b_st.site = &if_emlrtRSI;
  if ((1 <= x->size[1]) && (x->size[1] > 2147483646)) {
    c_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx - 1)) {
    if (x->data[ii]) {
      idx++;
      i->data[idx - 1] = ii + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  if (idx > x->size[1]) {
    emlrtErrorWithMessageIdR2018a(&st, &r_emlrtRTEI,
      "Coder:builtins:AssertionFailed", "Coder:builtins:AssertionFailed", 0);
  }

  if (x->size[1] == 1) {
    if (idx == 0) {
      i->size[0] = 1;
      i->size[1] = 0;
    }
  } else {
    ii = i->size[0] * i->size[1];
    if (1 > idx) {
      i->size[1] = 0;
    } else {
      i->size[1] = idx;
    }

    emxEnsureCapacity_int32_T(&st, i, ii, &od_emlrtRTEI);
  }
}

/* End of code generation (find.c) */
