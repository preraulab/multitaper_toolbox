/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * blockedSummation.c
 *
 * Code generation for function 'blockedSummation'
 *
 */

/* Include files */
#include "blockedSummation.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo ye_emlrtRSI = { 81, /* lineNo */
  "blockedSummation",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo af_emlrtRSI = { 142,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo bf_emlrtRSI = { 159,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo cf_emlrtRSI = { 161,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo df_emlrtRSI = { 173,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo ef_emlrtRSI = { 176,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo ff_emlrtRSI = { 194,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo gf_emlrtRSI = { 196,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo hf_emlrtRSI = { 207,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRTEInfo tc_emlrtRTEI = { 58,/* lineNo */
  9,                                   /* colNo */
  "blockedSummation",                  /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
};

static emlrtRTEInfo uc_emlrtRTEI = { 81,/* lineNo */
  13,                                  /* colNo */
  "blockedSummation",                  /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
};

static emlrtRTEInfo vc_emlrtRTEI = { 122,/* lineNo */
  1,                                   /* colNo */
  "blockedSummation",                  /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
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
  emxInit_real32_T(sp, &bsum, 1, &vc_emlrtRTEI, true);
  if ((x->size[0] == 0) || (x->size[1] == 0) || (vlen == 0)) {
    sz_idx_0 = (uint32_T)x->size[0];
    ib = y->size[0];
    y->size[0] = (int32_T)sz_idx_0;
    emxEnsureCapacity_real32_T(sp, y, ib, &tc_emlrtRTEI);
    firstBlockLength = (int32_T)sz_idx_0;
    for (ib = 0; ib < firstBlockLength; ib++) {
      y->data[ib] = 0.0F;
    }
  } else {
    st.site = &ye_emlrtRSI;
    vstride = x->size[0];
    bvstride = x->size[0] << 10;
    ib = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&st, y, ib, &uc_emlrtRTEI);
    ib = bsum->size[0];
    bsum->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&st, bsum, ib, &uc_emlrtRTEI);
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

    b_st.site = &af_emlrtRSI;
    if (x->size[0] > 2147483646) {
      c_st.site = &q_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (xj = 0; xj < vstride; xj++) {
      y->data[xj] = x->data[xj];
      bsum->data[xj] = 0.0F;
    }

    b_st.site = &bf_emlrtRSI;
    for (k = 2; k <= firstBlockLength; k++) {
      xoffset = (k - 1) * vstride;
      b_st.site = &cf_emlrtRSI;
      if (vstride > 2147483646) {
        c_st.site = &q_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      for (xj = 0; xj < vstride; xj++) {
        y->data[xj] += x->data[xoffset + xj];
      }
    }

    b_st.site = &df_emlrtRSI;
    for (ib = 2; ib <= nblocks; ib++) {
      firstBlockLength = (ib - 1) * bvstride;
      b_st.site = &ef_emlrtRSI;
      if (vstride > 2147483646) {
        c_st.site = &q_emlrtRSI;
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

      b_st.site = &ff_emlrtRSI;
      if ((2 <= hi) && (hi > 2147483646)) {
        c_st.site = &q_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      for (k = 2; k <= hi; k++) {
        xoffset = firstBlockLength + (k - 1) * vstride;
        b_st.site = &gf_emlrtRSI;
        for (xj = 0; xj < vstride; xj++) {
          bsum->data[xj] += x->data[xoffset + xj];
        }
      }

      b_st.site = &hf_emlrtRSI;
      for (xj = 0; xj < vstride; xj++) {
        y->data[xj] += bsum->data[xj];
      }
    }
  }

  emxFree_real32_T(&bsum);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (blockedSummation.c) */
