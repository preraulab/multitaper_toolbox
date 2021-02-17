/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * blockedSummation.cpp
 *
 * Code generation for function 'blockedSummation'
 *
 */

/* Include files */
#include "blockedSummation.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo qf_emlrtRSI = { 176,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo rf_emlrtRSI = { 196,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo sf_emlrtRSI = { 207,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRTEInfo fe_emlrtRTEI = { 58,/* lineNo */
  9,                                   /* colNo */
  "blockedSummation",                  /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
};

static emlrtRTEInfo ge_emlrtRTEI = { 122,/* lineNo */
  1,                                   /* colNo */
  "blockedSummation",                  /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
};

/* Function Definitions */
void blockedSummation(const emlrtStack *sp, const emxArray_real32_T *x, int32_T
                      vlen, emxArray_real32_T *y)
{
  emxArray_real32_T *bsum;
  uint32_T sz_idx_0;
  int32_T vstride;
  int32_T ib;
  int32_T bvstride;
  int32_T firstBlockLength;
  int32_T nblocks;
  int32_T lastBlockLength;
  int32_T xj;
  int32_T k;
  int32_T xoffset;
  int32_T hi;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real32_T(sp, &bsum, 1, &ge_emlrtRTEI, true);
  if ((x->size[1] == 0) || (vlen == 0)) {
    sz_idx_0 = static_cast<uint32_T>(x->size[0]);
    ib = y->size[0];
    y->size[0] = static_cast<int32_T>(sz_idx_0);
    emxEnsureCapacity_real32_T(sp, y, ib, &fe_emlrtRTEI);
    firstBlockLength = static_cast<int32_T>(sz_idx_0);
    for (ib = 0; ib < firstBlockLength; ib++) {
      y->data[ib] = 0.0F;
    }
  } else {
    st.site = &if_emlrtRSI;
    vstride = x->size[0];
    bvstride = x->size[0] << 10;
    ib = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&st, y, ib, &de_emlrtRTEI);
    ib = bsum->size[0];
    bsum->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&st, bsum, ib, &de_emlrtRTEI);
    if (vlen <= 1024) {
      firstBlockLength = vlen;
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = vlen / 1024;
      lastBlockLength = vlen - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }

    b_st.site = &jf_emlrtRSI;
    if ((1 <= x->size[0]) && (x->size[0] > 2147483646)) {
      c_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (xj = 0; xj < vstride; xj++) {
      y->data[xj] = x->data[xj];
      bsum->data[xj] = 0.0F;
    }

    b_st.site = &nf_emlrtRSI;
    for (k = 2; k <= firstBlockLength; k++) {
      xoffset = (k - 1) * vstride;
      b_st.site = &kf_emlrtRSI;
      if ((1 <= vstride) && (vstride > 2147483646)) {
        c_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      for (xj = 0; xj < vstride; xj++) {
        y->data[xj] += x->data[xoffset + xj];
      }
    }

    b_st.site = &of_emlrtRSI;
    for (ib = 2; ib <= nblocks; ib++) {
      firstBlockLength = (ib - 1) * bvstride;
      b_st.site = &qf_emlrtRSI;
      if ((1 <= vstride) && (vstride > 2147483646)) {
        c_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      for (xj = 0; xj < vstride; xj++) {
        bsum->data[xj] = x->data[firstBlockLength + xj];
      }

      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }

      b_st.site = &pf_emlrtRSI;
      if ((2 <= hi) && (hi > 2147483646)) {
        c_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      for (k = 2; k <= hi; k++) {
        xoffset = firstBlockLength + (k - 1) * vstride;
        b_st.site = &rf_emlrtRSI;
        for (xj = 0; xj < vstride; xj++) {
          bsum->data[xj] += x->data[xoffset + xj];
        }
      }

      b_st.site = &sf_emlrtRSI;
      for (xj = 0; xj < vstride; xj++) {
        y->data[xj] += bsum->data[xj];
      }
    }
  }

  emxFree_real32_T(&bsum);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (blockedSummation.cpp) */
