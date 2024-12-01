/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sumMatrixIncludeNaN.c
 *
 * Code generation for function 'sumMatrixIncludeNaN'
 *
 */

/* Include files */
#include "sumMatrixIncludeNaN.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo ee_emlrtRSI = {
    178,          /* lineNo */
    "sumColumnB", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo fe_emlrtRSI = {
    183,          /* lineNo */
    "sumColumnB", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo ge_emlrtRSI = {
    189,          /* lineNo */
    "sumColumnB", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo he_emlrtRSI = {
    210,         /* lineNo */
    "sumColumn", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

/* Function Definitions */
real32_T b_sumColumnB(const emlrtStack *sp, const emxArray_real32_T *x,
                      int32_T col, int32_T vlen, int32_T vstart)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  int32_T b_k;
  int32_T k;
  const real32_T *x_data;
  real32_T y;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  x_data = x->data;
  if (vlen <= 1024) {
    int32_T i0;
    st.site = &ee_emlrtRSI;
    i0 = vstart + (col - 1) * x->size[0];
    y = x_data[i0 - 1];
    b_st.site = &he_emlrtRSI;
    if (vlen - 1 > 2147483646) {
      c_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (k = 0; k <= vlen - 2; k++) {
      y += x_data[i0 + k];
    }
  } else {
    int32_T i0;
    int32_T i0_tmp;
    int32_T inb;
    int32_T nfb;
    real32_T b_y;
    nfb = vlen / 1024;
    inb = nfb << 10;
    i0_tmp = (col - 1) * x->size[0];
    i0 = vstart + i0_tmp;
    y = x_data[i0 - 1];
    for (k = 0; k < 1023; k++) {
      y += x_data[i0 + k];
    }
    st.site = &fe_emlrtRSI;
    for (k = 2; k <= nfb; k++) {
      i0 = (vstart + ((k - 1) << 10)) + i0_tmp;
      b_y = x_data[i0 - 1];
      for (b_k = 0; b_k < 1023; b_k++) {
        b_y += x_data[i0 + b_k];
      }
      y += b_y;
    }
    if (vlen > inb) {
      st.site = &ge_emlrtRSI;
      i0 = (vstart + inb) + i0_tmp;
      b_y = x_data[i0 - 1];
      nfb = vlen - inb;
      b_st.site = &he_emlrtRSI;
      for (k = 0; k <= nfb - 2; k++) {
        b_y += x_data[i0 + k];
      }
      y += b_y;
    }
  }
  return y;
}

real32_T sumColumnB(const emlrtStack *sp, const emxArray_real32_T *x,
                    int32_T col, int32_T vlen)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  int32_T b_k;
  int32_T k;
  const real32_T *x_data;
  real32_T y;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  x_data = x->data;
  if (vlen <= 1024) {
    int32_T i0;
    int32_T nfb;
    st.site = &ee_emlrtRSI;
    i0 = (col - 1) * x->size[0];
    y = x_data[i0];
    b_st.site = &he_emlrtRSI;
    if (vlen - 1 > 2147483646) {
      c_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    nfb = (uint16_T)(vlen - 1);
    for (k = 0; k < nfb; k++) {
      y += x_data[(i0 + k) + 1];
    }
  } else {
    int32_T i0;
    int32_T i0_tmp;
    int32_T inb;
    int32_T nfb;
    real32_T b_y;
    nfb = vlen / 1024;
    inb = nfb << 10;
    i0_tmp = (col - 1) * x->size[0];
    y = x_data[i0_tmp];
    for (k = 0; k < 1023; k++) {
      y += x_data[(i0_tmp + k) + 1];
    }
    st.site = &fe_emlrtRSI;
    for (k = 2; k <= nfb; k++) {
      i0 = ((k - 1) << 10) + i0_tmp;
      b_y = x_data[i0];
      for (b_k = 0; b_k < 1023; b_k++) {
        b_y += x_data[(i0 + b_k) + 1];
      }
      y += b_y;
    }
    if (vlen > inb) {
      st.site = &ge_emlrtRSI;
      i0 = (inb + i0_tmp) + 1;
      b_y = x_data[i0 - 1];
      nfb = (vlen - inb) - 2;
      b_st.site = &he_emlrtRSI;
      for (k = 0; k <= nfb; k++) {
        b_y += x_data[i0 + k];
      }
      y += b_y;
    }
  }
  return y;
}

real32_T sumColumnB4(const emxArray_real32_T *x, int32_T col, int32_T vstart)
{
  int32_T i1;
  int32_T k;
  const real32_T *x_data;
  real32_T psum2;
  real32_T psum3;
  real32_T psum4;
  real32_T y;
  x_data = x->data;
  i1 = vstart + (col - 1) * x->size[0];
  y = x_data[i1 - 1];
  psum2 = x_data[i1 + 1023];
  psum3 = x_data[i1 + 2047];
  psum4 = x_data[i1 + 3071];
  for (k = 0; k < 1023; k++) {
    int32_T psum1_tmp;
    psum1_tmp = i1 + k;
    y += x_data[psum1_tmp];
    psum2 += x_data[psum1_tmp + 1024];
    psum3 += x_data[psum1_tmp + 2048];
    psum4 += x_data[psum1_tmp + 3072];
  }
  return (y + psum2) + (psum3 + psum4);
}

/* End of code generation (sumMatrixIncludeNaN.c) */
