/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.cpp
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo mf_emlrtRSI = { 137,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRTEInfo c_emlrtRTEI = { 43,/* lineNo */
  23,                                  /* colNo */
  "sumprod",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pName */
};

static emlrtRTEInfo ee_emlrtRTEI = { 20,/* lineNo */
  1,                                   /* colNo */
  "sum",                               /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/sum.m"/* pName */
};

/* Function Definitions */
void sum(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T *y)
{
  uint32_T sz_idx_1;
  int32_T npages;
  int32_T nblocks;
  int32_T firstBlockLength;
  int32_T lastBlockLength;
  int32_T xi;
  int32_T xpageoffset;
  int32_T k;
  int32_T ib;
  int32_T xblockoffset;
  real32_T bsum;
  int32_T hi;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &lf_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  if ((x->size[0] == 1) && (x->size[1] != 1)) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:toolbox:autoDimIncompatibility",
      "Coder:toolbox:autoDimIncompatibility", 0);
  }

  b_st.site = &jb_emlrtRSI;
  c_st.site = &hf_emlrtRSI;
  if (x->size[0] == 0) {
    sz_idx_1 = static_cast<uint32_T>(x->size[1]);
    nblocks = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = static_cast<int32_T>(sz_idx_1);
    emxEnsureCapacity_real32_T(&c_st, y, nblocks, &ee_emlrtRTEI);
    firstBlockLength = static_cast<int32_T>(sz_idx_1);
    for (nblocks = 0; nblocks < firstBlockLength; nblocks++) {
      y->data[nblocks] = 0.0F;
    }
  } else {
    d_st.site = &if_emlrtRSI;
    npages = x->size[1];
    nblocks = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real32_T(&d_st, y, nblocks, &de_emlrtRTEI);
    if (x->size[0] <= 1024) {
      firstBlockLength = x->size[0];
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = x->size[0] / 1024;
      lastBlockLength = x->size[0] - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }

    e_st.site = &mf_emlrtRSI;
    if ((1 <= x->size[1]) && (x->size[1] > 2147483646)) {
      f_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&f_st);
    }

    for (xi = 0; xi < npages; xi++) {
      xpageoffset = xi * x->size[0];
      y->data[xi] = x->data[xpageoffset];
      e_st.site = &nf_emlrtRSI;
      for (k = 2; k <= firstBlockLength; k++) {
        y->data[xi] += x->data[(xpageoffset + k) - 1];
      }

      e_st.site = &of_emlrtRSI;
      for (ib = 2; ib <= nblocks; ib++) {
        xblockoffset = xpageoffset + ((ib - 1) << 10);
        bsum = x->data[xblockoffset];
        if (ib == nblocks) {
          hi = lastBlockLength;
        } else {
          hi = 1024;
        }

        e_st.site = &pf_emlrtRSI;
        if ((2 <= hi) && (hi > 2147483646)) {
          f_st.site = &ab_emlrtRSI;
          check_forloop_overflow_error(&f_st);
        }

        for (k = 2; k <= hi; k++) {
          bsum += x->data[(xblockoffset + k) - 1];
        }

        y->data[xi] += bsum;
      }
    }
  }
}

/* End of code generation (sum.cpp) */
