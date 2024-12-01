/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * colon.c
 *
 * Code generation for function 'colon'
 *
 */

/* Include files */
#include "colon.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo x_emlrtRSI = {
    319,               /* lineNo */
    "eml_float_colon", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo pb_emlrtRSI = {
    148,                            /* lineNo */
    "eml_integer_colon_dispatcher", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo qb_emlrtRSI = {
    171,                        /* lineNo */
    "eml_signed_integer_colon", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo rb_emlrtRSI = {
    176,                        /* lineNo */
    "eml_signed_integer_colon", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRTEInfo d_emlrtRTEI = {
    419,               /* lineNo */
    15,                /* colNo */
    "assert_pmaxsize", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pName
                                                                          */
};

static emlrtRTEInfo wc_emlrtRTEI = {
    320,     /* lineNo */
    20,      /* colNo */
    "colon", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pName
                                                                          */
};

static emlrtRTEInfo bd_emlrtRTEI = {
    172,     /* lineNo */
    20,      /* colNo */
    "colon", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pName
                                                                          */
};

/* Function Definitions */
void b_eml_float_colon(const emlrtStack *sp, real_T d, real_T b,
                       emxArray_real_T *y)
{
  emlrtStack st;
  real_T apnd;
  real_T cdiff;
  real_T ndbl;
  real_T *y_data;
  int32_T k;
  int32_T n;
  int32_T nm1d2;
  st.prev = sp;
  st.tls = sp->tls;
  ndbl = muDoubleScalarFloor(b / d + 0.5);
  apnd = ndbl * d;
  if (d > 0.0) {
    cdiff = apnd - b;
  } else {
    cdiff = b - apnd;
  }
  if (muDoubleScalarAbs(cdiff) <
      4.4408920985006262E-16 * muDoubleScalarMax(0.0, muDoubleScalarAbs(b))) {
    ndbl++;
    apnd = b;
  } else if (cdiff > 0.0) {
    apnd = (ndbl - 1.0) * d;
  } else {
    ndbl++;
  }
  if (ndbl >= 0.0) {
    n = (int32_T)ndbl;
  } else {
    n = 0;
  }
  st.site = &x_emlrtRSI;
  if (ndbl > 2.147483647E+9) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI, "Coder:MATLAB:pmaxsize",
                                  "Coder:MATLAB:pmaxsize", 0);
  }
  nm1d2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity_real_T(sp, y, nm1d2, &wc_emlrtRTEI);
  y_data = y->data;
  if (n > 0) {
    y_data[0] = 0.0;
    if (n > 1) {
      y_data[n - 1] = apnd;
      nm1d2 = (n - 1) / 2;
      for (k = 0; k <= nm1d2 - 2; k++) {
        ndbl = ((real_T)k + 1.0) * d;
        y_data[k + 1] = ndbl;
        y_data[(n - k) - 2] = apnd - ndbl;
      }
      if (nm1d2 << 1 == n - 1) {
        y_data[nm1d2] = apnd / 2.0;
      } else {
        ndbl = (real_T)nm1d2 * d;
        y_data[nm1d2] = ndbl;
        y_data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }
}

void eml_float_colon(const emlrtStack *sp, real_T d, real_T b,
                     emxArray_real_T *y)
{
  emlrtStack st;
  real_T apnd;
  real_T cdiff;
  real_T ndbl;
  real_T *y_data;
  int32_T k;
  int32_T n;
  int32_T nm1d2;
  st.prev = sp;
  st.tls = sp->tls;
  ndbl = muDoubleScalarFloor((b - 1.0) / d + 0.5);
  apnd = ndbl * d + 1.0;
  if (d > 0.0) {
    cdiff = apnd - b;
  } else {
    cdiff = b - apnd;
  }
  if (muDoubleScalarAbs(cdiff) <
      4.4408920985006262E-16 * muDoubleScalarMax(1.0, muDoubleScalarAbs(b))) {
    ndbl++;
    apnd = b;
  } else if (cdiff > 0.0) {
    apnd = (ndbl - 1.0) * d + 1.0;
  } else {
    ndbl++;
  }
  if (ndbl >= 0.0) {
    n = (int32_T)ndbl;
  } else {
    n = 0;
  }
  st.site = &x_emlrtRSI;
  if (ndbl > 2.147483647E+9) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI, "Coder:MATLAB:pmaxsize",
                                  "Coder:MATLAB:pmaxsize", 0);
  }
  nm1d2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity_real_T(sp, y, nm1d2, &wc_emlrtRTEI);
  y_data = y->data;
  if (n > 0) {
    y_data[0] = 1.0;
    if (n > 1) {
      y_data[n - 1] = apnd;
      nm1d2 = (n - 1) / 2;
      for (k = 0; k <= nm1d2 - 2; k++) {
        ndbl = ((real_T)k + 1.0) * d;
        y_data[k + 1] = ndbl + 1.0;
        y_data[(n - k) - 2] = apnd - ndbl;
      }
      if (nm1d2 << 1 == n - 1) {
        y_data[nm1d2] = (apnd + 1.0) / 2.0;
      } else {
        ndbl = (real_T)nm1d2 * d;
        y_data[nm1d2] = ndbl + 1.0;
        y_data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }
}

void eml_integer_colon_dispatcher(const emlrtStack *sp, int32_T b,
                                  emxArray_int32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  int32_T k;
  int32_T n;
  int32_T yk;
  int32_T *y_data;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &pb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  b_st.site = &qb_emlrtRSI;
  if (b < 0) {
    n = 0;
  } else {
    n = b + 1;
  }
  yk = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity_int32_T(&st, y, yk, &bd_emlrtRTEI);
  y_data = y->data;
  if (n > 0) {
    y_data[0] = 0;
    yk = 0;
    b_st.site = &rb_emlrtRSI;
    if (n > 2147483646) {
      c_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (k = 2; k <= n; k++) {
      yk++;
      y_data[k - 1] = yk;
    }
  }
}

/* End of code generation (colon.c) */
