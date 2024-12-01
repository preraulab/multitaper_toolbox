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
static emlrtRSInfo pe_emlrtRSI = {
    138,        /* lineNo */
    "eml_find", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                           */
};

static emlrtRSInfo qe_emlrtRSI = {
    376,                  /* lineNo */
    "find_first_indices", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                           */
};

static emlrtRTEInfo p_emlrtRTEI = {
    386,                  /* lineNo */
    1,                    /* colNo */
    "find_first_indices", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/find.m" /* pName
                                                                           */
};

static emlrtRTEInfo kd_emlrtRTEI = {
    358,    /* lineNo */
    24,     /* colNo */
    "find", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/find.m" /* pName
                                                                           */
};

static emlrtRTEInfo ld_emlrtRTEI = {
    138,    /* lineNo */
    9,      /* colNo */
    "find", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/find.m" /* pName
                                                                           */
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
  int32_T nx_tmp;
  int32_T *i_data;
  const boolean_T *x_data;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  x_data = x->data;
  nx_tmp = x->size[1];
  st.site = &pe_emlrtRSI;
  idx = 0;
  ii = i->size[0] * i->size[1];
  i->size[0] = 1;
  i->size[1] = x->size[1];
  emxEnsureCapacity_int32_T(&st, i, ii, &kd_emlrtRTEI);
  i_data = i->data;
  b_st.site = &qe_emlrtRSI;
  if (x->size[1] > 2147483646) {
    c_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx_tmp - 1)) {
    if (x_data[ii]) {
      idx++;
      i_data[idx - 1] = ii + 1;
      if (idx >= nx_tmp) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }
  if (idx > x->size[1]) {
    emlrtErrorWithMessageIdR2018a(&st, &p_emlrtRTEI,
                                  "Coder:builtins:AssertionFailed",
                                  "Coder:builtins:AssertionFailed", 0);
  }
  if (x->size[1] == 1) {
    if (idx == 0) {
      i->size[0] = 1;
      i->size[1] = 0;
    }
  } else {
    ii = i->size[0] * i->size[1];
    if (idx < 1) {
      i->size[1] = 0;
    } else {
      i->size[1] = idx;
    }
    emxEnsureCapacity_int32_T(&st, i, ii, &ld_emlrtRTEI);
  }
}

/* End of code generation (find.c) */
