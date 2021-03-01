/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo xe_emlrtRSI = { 20, /* lineNo */
  "sum",                               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

static emlrtRSInfo ye_emlrtRSI = { 144,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo af_emlrtRSI = { 166,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo bf_emlrtRSI = { 180,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo cf_emlrtRSI = { 201,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo df_emlrtRSI = { 183,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo ef_emlrtRSI = { 203,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRSInfo ff_emlrtRSI = { 214,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

static emlrtRTEInfo c_emlrtRTEI = { 46,/* lineNo */
  23,                                  /* colNo */
  "sumprod",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pName */
};

static emlrtRTEInfo kd_emlrtRTEI = { 20,/* lineNo */
  1,                                   /* colNo */
  "sum",                               /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/sum.m"/* pName */
};

static emlrtRTEInfo ld_emlrtRTEI = { 129,/* lineNo */
  23,                                  /* colNo */
  "blockedSummation",                  /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
};

static emlrtRTEInfo md_emlrtRTEI = { 129,/* lineNo */
  1,                                   /* colNo */
  "blockedSummation",                  /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
};

/* Function Definitions */
void b_sum(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T
           *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_real32_T *bsum;
  int32_T bvstride;
  int32_T firstBlockLength;
  int32_T hi;
  int32_T ib;
  int32_T k;
  int32_T lastBlockLength;
  int32_T nblocks;
  int32_T vstride;
  int32_T xj;
  int32_T xoffset;
  st.prev = sp;
  st.tls = sp->tls;
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
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  st.site = &xe_emlrtRSI;
  b_st.site = &jb_emlrtRSI;
  c_st.site = &te_emlrtRSI;
  if (x->size[1] == 0) {
    hi = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&c_st, y, hi, &kd_emlrtRTEI);
    firstBlockLength = x->size[0];
    for (hi = 0; hi < firstBlockLength; hi++) {
      y->data[hi] = 0.0F;
    }
  } else {
    emxInit_real32_T(&c_st, &bsum, 1, &md_emlrtRTEI, true);
    d_st.site = &ue_emlrtRSI;
    vstride = x->size[0];
    bvstride = x->size[0] << 10;
    hi = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&d_st, y, hi, &jd_emlrtRTEI);
    hi = bsum->size[0];
    bsum->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&d_st, bsum, hi, &ld_emlrtRTEI);
    if (x->size[1] <= 1024) {
      firstBlockLength = x->size[1];
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = x->size[1] / 1024;
      lastBlockLength = x->size[1] - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }

    e_st.site = &ve_emlrtRSI;
    if ((1 <= x->size[0]) && (x->size[0] > 2147483646)) {
      f_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&f_st);
    }

    for (xj = 0; xj < vstride; xj++) {
      y->data[xj] = x->data[xj];
      bsum->data[xj] = 0.0F;
    }

    e_st.site = &af_emlrtRSI;
    for (k = 2; k <= firstBlockLength; k++) {
      xoffset = (k - 1) * vstride;
      e_st.site = &we_emlrtRSI;
      if ((1 <= vstride) && (vstride > 2147483646)) {
        f_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      for (xj = 0; xj < vstride; xj++) {
        y->data[xj] += x->data[xoffset + xj];
      }
    }

    e_st.site = &bf_emlrtRSI;
    for (ib = 2; ib <= nblocks; ib++) {
      firstBlockLength = (ib - 1) * bvstride;
      e_st.site = &df_emlrtRSI;
      if ((1 <= vstride) && (vstride > 2147483646)) {
        f_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      for (xj = 0; xj < vstride; xj++) {
        bsum->data[xj] = x->data[firstBlockLength + xj];
      }

      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }

      e_st.site = &cf_emlrtRSI;
      if ((2 <= hi) && (hi > 2147483646)) {
        f_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      for (k = 2; k <= hi; k++) {
        xoffset = firstBlockLength + (k - 1) * vstride;
        e_st.site = &ef_emlrtRSI;
        for (xj = 0; xj < vstride; xj++) {
          bsum->data[xj] += x->data[xoffset + xj];
        }
      }

      e_st.site = &ff_emlrtRSI;
      for (xj = 0; xj < vstride; xj++) {
        y->data[xj] += bsum->data[xj];
      }
    }

    emxFree_real32_T(&bsum);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void sum(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  int32_T firstBlockLength;
  int32_T hi;
  int32_T ib;
  int32_T k;
  int32_T lastBlockLength;
  int32_T nblocks;
  int32_T npages;
  int32_T xblockoffset;
  int32_T xi;
  int32_T xpageoffset;
  real32_T bsum;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &xe_emlrtRSI;
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
  if (((x->size[0] != 1) || (x->size[1] != 1)) && (x->size[0] == 1)) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:toolbox:autoDimIncompatibility",
      "Coder:toolbox:autoDimIncompatibility", 0);
  }

  b_st.site = &jb_emlrtRSI;
  c_st.site = &te_emlrtRSI;
  if (x->size[0] == 0) {
    nblocks = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real32_T(&c_st, y, nblocks, &kd_emlrtRTEI);
    firstBlockLength = x->size[1];
    for (nblocks = 0; nblocks < firstBlockLength; nblocks++) {
      y->data[nblocks] = 0.0F;
    }
  } else {
    d_st.site = &ue_emlrtRSI;
    npages = x->size[1];
    nblocks = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real32_T(&d_st, y, nblocks, &jd_emlrtRTEI);
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

    e_st.site = &ye_emlrtRSI;
    if ((1 <= x->size[1]) && (x->size[1] > 2147483646)) {
      f_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&f_st);
    }

    for (xi = 0; xi < npages; xi++) {
      xpageoffset = xi * x->size[0];
      y->data[xi] = x->data[xpageoffset];
      e_st.site = &af_emlrtRSI;
      for (k = 2; k <= firstBlockLength; k++) {
        y->data[xi] += x->data[(xpageoffset + k) - 1];
      }

      e_st.site = &bf_emlrtRSI;
      for (ib = 2; ib <= nblocks; ib++) {
        xblockoffset = xpageoffset + ((ib - 1) << 10);
        bsum = x->data[xblockoffset];
        if (ib == nblocks) {
          hi = lastBlockLength;
        } else {
          hi = 1024;
        }

        e_st.site = &cf_emlrtRSI;
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

/* End of code generation (sum.c) */
