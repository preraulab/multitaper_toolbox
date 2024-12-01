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
#include "sumMatrixIncludeNaN.h"

/* Variable Definitions */
static emlrtRSInfo vd_emlrtRSI = {
    20,    /* lineNo */
    "sum", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/sum.m" /* pathName
                                                                            */
};

static emlrtRSInfo wd_emlrtRSI = {
    107,                /* lineNo */
    "blockedSummation", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo xd_emlrtRSI = {
    22,                    /* lineNo */
    "sumMatrixIncludeNaN", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo yd_emlrtRSI = {
    41,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo ae_emlrtRSI = {
    42,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo be_emlrtRSI = {
    50,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo ce_emlrtRSI = {
    53,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo de_emlrtRSI = {
    57,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo ie_emlrtRSI = {
    190,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo je_emlrtRSI = {
    204,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo ke_emlrtRSI = {
    207,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo le_emlrtRSI = {
    225,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo me_emlrtRSI = {
    227,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo ne_emlrtRSI = {
    238,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRTEInfo c_emlrtRTEI = {
    46,        /* lineNo */
    23,        /* colNo */
    "sumprod", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumprod.m" /* pName */
};

static emlrtRTEInfo o_emlrtRTEI = {
    76,        /* lineNo */
    9,         /* colNo */
    "sumprod", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumprod.m" /* pName */
};

static emlrtRTEInfo gd_emlrtRTEI = {
    20,    /* lineNo */
    1,     /* colNo */
    "sum", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/sum.m" /* pName
                                                                            */
};

static emlrtRTEInfo hd_emlrtRTEI = {
    35,                    /* lineNo */
    20,                    /* colNo */
    "sumMatrixIncludeNaN", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pName */
};

static emlrtRTEInfo id_emlrtRTEI = {
    153,                /* lineNo */
    23,                 /* colNo */
    "blockedSummation", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pName */
};

static emlrtRTEInfo jd_emlrtRTEI = {
    153,                /* lineNo */
    1,                  /* colNo */
    "blockedSummation", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pName */
};

/* Function Definitions */
void b_sum(const emlrtStack *sp, const emxArray_real32_T *x,
           emxArray_real32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_real32_T *bsum;
  int32_T hi;
  int32_T ib;
  int32_T k;
  int32_T xj;
  const real32_T *x_data;
  real32_T *bsum_data;
  real32_T *y_data;
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
  x_data = x->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  st.site = &vd_emlrtRSI;
  b_st.site = &eb_emlrtRSI;
  c_st.site = &pd_emlrtRSI;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    int32_T firstBlockLength;
    firstBlockLength = x->size[0];
    hi = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&c_st, y, hi, &gd_emlrtRTEI);
    y_data = y->data;
    for (hi = 0; hi < firstBlockLength; hi++) {
      y_data[hi] = 0.0F;
    }
  } else {
    int32_T bvstride;
    int32_T firstBlockLength;
    int32_T lastBlockLength;
    int32_T nblocks;
    int32_T vstride_tmp;
    int32_T xoffset;
    d_st.site = &qd_emlrtRSI;
    vstride_tmp = x->size[0];
    bvstride = x->size[0] << 10;
    hi = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&d_st, y, hi, &fd_emlrtRTEI);
    y_data = y->data;
    emxInit_real32_T(&d_st, &bsum, 1, &jd_emlrtRTEI);
    hi = bsum->size[0];
    bsum->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&d_st, bsum, hi, &id_emlrtRTEI);
    bsum_data = bsum->data;
    if (x->size[1] <= 1024) {
      firstBlockLength = x->size[1];
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = (int32_T)((uint32_T)x->size[1] >> 10);
      lastBlockLength = x->size[1] - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    e_st.site = &rd_emlrtRSI;
    if (x->size[0] > 2147483646) {
      f_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&f_st);
    }
    for (xj = 0; xj < vstride_tmp; xj++) {
      y_data[xj] = x_data[xj];
      bsum_data[xj] = 0.0F;
    }
    e_st.site = &ie_emlrtRSI;
    for (k = 2; k <= firstBlockLength; k++) {
      xoffset = (k - 1) * vstride_tmp;
      e_st.site = &sd_emlrtRSI;
      if (vstride_tmp > 2147483646) {
        f_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }
      for (xj = 0; xj < vstride_tmp; xj++) {
        y_data[xj] += x_data[xoffset + xj];
      }
    }
    e_st.site = &je_emlrtRSI;
    for (ib = 2; ib <= nblocks; ib++) {
      firstBlockLength = (ib - 1) * bvstride;
      e_st.site = &ke_emlrtRSI;
      if (vstride_tmp > 2147483646) {
        f_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }
      for (xj = 0; xj < vstride_tmp; xj++) {
        bsum_data[xj] = x_data[firstBlockLength + xj];
      }
      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      e_st.site = &le_emlrtRSI;
      for (k = 2; k <= hi; k++) {
        xoffset = firstBlockLength + (k - 1) * vstride_tmp;
        e_st.site = &me_emlrtRSI;
        for (xj = 0; xj < vstride_tmp; xj++) {
          bsum_data[xj] += x_data[xoffset + xj];
        }
      }
      e_st.site = &ne_emlrtRSI;
      for (xj = 0; xj < vstride_tmp; xj++) {
        y_data[xj] += bsum_data[xj];
      }
    }
    emxFree_real32_T(&d_st, &bsum);
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

void sum(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  int32_T col;
  int32_T ib;
  int32_T inb;
  real32_T *y_data;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &vd_emlrtRSI;
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
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  if (((x->size[0] != 1) || (x->size[1] != 1)) && (x->size[0] == 1)) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
                                  "Coder:toolbox:autoDimIncompatibility",
                                  "Coder:toolbox:autoDimIncompatibility", 0);
  }
  if ((x->size[0] == 0) && (x->size[1] == 0)) {
    emlrtErrorWithMessageIdR2018a(&st, &o_emlrtRTEI,
                                  "Coder:toolbox:UnsupportedSpecialEmpty",
                                  "Coder:toolbox:UnsupportedSpecialEmpty", 0);
  }
  b_st.site = &eb_emlrtRSI;
  c_st.site = &pd_emlrtRSI;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    int32_T nfb;
    inb = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real32_T(&c_st, y, inb, &gd_emlrtRTEI);
    y_data = y->data;
    nfb = x->size[1];
    for (inb = 0; inb < nfb; inb++) {
      y_data[inb] = 0.0F;
    }
  } else {
    int32_T i;
    d_st.site = &wd_emlrtRSI;
    e_st.site = &xd_emlrtRSI;
    inb = y->size[0] * y->size[1];
    y->size[0] = 1;
    i = x->size[1];
    y->size[1] = x->size[1];
    emxEnsureCapacity_real32_T(&e_st, y, inb, &hd_emlrtRTEI);
    y_data = y->data;
    if (x->size[0] < 4096) {
      f_st.site = &yd_emlrtRSI;
      if (x->size[1] > 2147483646) {
        g_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&g_st);
      }
      for (col = 0; col < i; col++) {
        f_st.site = &ae_emlrtRSI;
        y_data[col] = sumColumnB(&f_st, x, col + 1, x->size[0]);
      }
    } else {
      int32_T nfb;
      int32_T nleft;
      nfb = (int32_T)((uint32_T)x->size[0] >> 12);
      inb = nfb << 12;
      nleft = x->size[0] - inb;
      f_st.site = &be_emlrtRSI;
      if (x->size[1] > 2147483646) {
        g_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&g_st);
      }
      for (col = 0; col < i; col++) {
        real32_T s;
        s = sumColumnB4(x, col + 1, 1);
        f_st.site = &ce_emlrtRSI;
        for (ib = 2; ib <= nfb; ib++) {
          s += sumColumnB4(x, col + 1, ((ib - 1) << 12) + 1);
        }
        if (nleft > 0) {
          f_st.site = &de_emlrtRSI;
          s += b_sumColumnB(&f_st, x, col + 1, nleft, inb + 1);
        }
        y_data[col] = s;
      }
    }
  }
}

/* End of code generation (sum.c) */
