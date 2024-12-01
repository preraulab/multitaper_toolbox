/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * detrend.c
 *
 * Code generation for function 'detrend'
 *
 */

/* Include files */
#include "detrend.h"
#include "colon.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "blas.h"
#include "lapacke.h"
#include "mwmathutil.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo lb_emlrtRSI =
    {
        102,       /* lineNo */
        "detrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI =
    {
        198,       /* lineNo */
        "detrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo nb_emlrtRSI = {
    28,      /* lineNo */
    "colon", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo ob_emlrtRSI = {
    82,      /* lineNo */
    "colon", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo ub_emlrtRSI =
    {
        208,                   /* lineNo */
        "chooseDetrendMethod", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo vb_emlrtRSI =
    {
        369,         /* lineNo */
        "CCDetrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo wb_emlrtRSI = {
    66,                    /* lineNo */
    "applyVectorFunction", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
    "applyVectorFunction.m" /* pathName */
};

static emlrtRSInfo xb_emlrtRSI = {
    63,                               /* lineNo */
    "function_handle/parenReference", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
    "function_handle.m" /* pathName */
};

static emlrtRSInfo yb_emlrtRSI =
    {
        368,                        /* lineNo */
        "@(x)subsum(x,ONE,endSeg)", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo ac_emlrtRSI =
    {
        268,      /* lineNo */
        "subsum", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo bc_emlrtRSI =
    {
        271,      /* lineNo */
        "subsum", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo cc_emlrtRSI =
    {
        278,      /* lineNo */
        "subsum", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo dc_emlrtRSI =
    {
        288,     /* lineNo */
        "sumab", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo ec_emlrtRSI =
    {
        210,                   /* lineNo */
        "chooseDetrendMethod", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo fc_emlrtRSI =
    {
        391,          /* lineNo */
        "CPPDetrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo gc_emlrtRSI =
    {
        401,          /* lineNo */
        "CPPDetrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo hc_emlrtRSI =
    {
        402,          /* lineNo */
        "CPPDetrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo ic_emlrtRSI =
    {
        406,          /* lineNo */
        "CPPDetrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo jc_emlrtRSI =
    {
        411,          /* lineNo */
        "CPPDetrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo kc_emlrtRSI =
    {
        414,          /* lineNo */
        "CPPDetrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo lc_emlrtRSI =
    {
        416,          /* lineNo */
        "CPPDetrend", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pathName */
};

static emlrtRSInfo mc_emlrtRSI =
    {
        61,        /* lineNo */
        "qrsolve", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
        "qrsolve.m" /* pathName */
};

static emlrtRSInfo nc_emlrtRSI =
    {
        85,        /* lineNo */
        "qrsolve", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
        "qrsolve.m" /* pathName */
};

static emlrtRSInfo oc_emlrtRSI = {
    63,       /* lineNo */
    "xgeqp3", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgeqp3.m" /* pathName */
};

static emlrtRSInfo pc_emlrtRSI = {
    138,            /* lineNo */
    "ceval_xgeqp3", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgeqp3.m" /* pathName */
};

static emlrtRSInfo qc_emlrtRSI = {
    143,            /* lineNo */
    "ceval_xgeqp3", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgeqp3.m" /* pathName */
};

static emlrtRSInfo rc_emlrtRSI = {
    148,            /* lineNo */
    "ceval_xgeqp3", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgeqp3.m" /* pathName */
};

static emlrtRSInfo sc_emlrtRSI = {
    151,            /* lineNo */
    "ceval_xgeqp3", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgeqp3.m" /* pathName */
};

static emlrtRSInfo tc_emlrtRSI =
    {
        119,         /* lineNo */
        "LSQFromQR", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
        "qrsolve.m" /* pathName */
};

static emlrtRSInfo uc_emlrtRSI =
    {
        128,         /* lineNo */
        "LSQFromQR", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
        "qrsolve.m" /* pathName */
};

static emlrtRSInfo vc_emlrtRSI =
    {
        138,         /* lineNo */
        "LSQFromQR", /* fcnName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
        "qrsolve.m" /* pathName */
};

static emlrtRSInfo wc_emlrtRSI = {
    40,         /* lineNo */
    "xunormqr", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xunormqr.m" /* pathName */
};

static emlrtRSInfo xc_emlrtRSI = {
    106,              /* lineNo */
    "ceval_xunormqr", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xunormqr.m" /* pathName */
};

static emlrtRSInfo yc_emlrtRSI = {
    94,                  /* lineNo */
    "eml_mtimes_helper", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/"
    "eml_mtimes_helper.m" /* pathName */
};

static emlrtRSInfo ad_emlrtRSI = {
    142,      /* lineNo */
    "mtimes", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "mtimes.m" /* pathName */
};

static emlrtRSInfo bd_emlrtRSI = {
    178,           /* lineNo */
    "mtimes_blas", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "mtimes.m" /* pathName */
};

static emlrtRTEInfo v_emlrtRTEI = {
    48,          /* lineNo */
    13,          /* colNo */
    "infocheck", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "infocheck.m" /* pName */
};

static emlrtRTEInfo w_emlrtRTEI = {
    45,          /* lineNo */
    13,          /* colNo */
    "infocheck", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "infocheck.m" /* pName */
};

static emlrtRTEInfo xc_emlrtRTEI =
    {
        141,       /* lineNo */
        5,         /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo
    yc_emlrtRTEI =
        {
            60,       /* lineNo */
            20,       /* colNo */
            "bsxfun", /* fName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/"
            "bsxfun.m" /* pName */
};

static emlrtRTEInfo ad_emlrtRTEI =
    {
        102,       /* lineNo */
        78,        /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo rd_emlrtRTEI =
    {
        102,       /* lineNo */
        71,        /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo sd_emlrtRTEI =
    {
        389,       /* lineNo */
        20,        /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo td_emlrtRTEI =
    {
        397,       /* lineNo */
        20,        /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo ud_emlrtRTEI = {
    1,        /* lineNo */
    32,       /* colNo */
    "xgeqp3", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgeqp3.m" /* pName */
};

static emlrtRTEInfo vd_emlrtRTEI =
    {
        85,        /* lineNo */
        26,        /* colNo */
        "qrsolve", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
        "qrsolve.m" /* pName */
};

static emlrtRTEInfo wd_emlrtRTEI =
    {
        119,       /* lineNo */
        5,         /* colNo */
        "qrsolve", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
        "qrsolve.m" /* pName */
};

static emlrtRTEInfo xd_emlrtRTEI = {
    218,      /* lineNo */
    20,       /* colNo */
    "mtimes", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "mtimes.m" /* pName */
};

static emlrtRTEInfo yd_emlrtRTEI =
    {
        102,       /* lineNo */
        1,         /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo ae_emlrtRTEI =
    {
        389,       /* lineNo */
        1,         /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo be_emlrtRTEI =
    {
        397,       /* lineNo */
        1,         /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo ce_emlrtRTEI =
    {
        1,         /* lineNo */
        14,        /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

static emlrtRTEInfo fe_emlrtRTEI =
    {
        416,       /* lineNo */
        5,         /* colNo */
        "detrend", /* fName */
        "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/"
        "detrend.m" /* pName */
};

/* Function Declarations */
static void minus(const emlrtStack *sp, emxArray_real32_T *in1,
                  const emxArray_real32_T *in2);

/* Function Definitions */
static void minus(const emlrtStack *sp, emxArray_real32_T *in1,
                  const emxArray_real32_T *in2)
{
  emxArray_real32_T *b_in1;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_1_0;
  const real32_T *in2_data;
  real32_T *b_in1_data;
  real32_T *in1_data;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real32_T(sp, &b_in1, 1, &fe_emlrtRTEI);
  if (in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = in2->size[0];
  }
  i = b_in1->size[0];
  b_in1->size[0] = loop_ub;
  emxEnsureCapacity_real32_T(sp, b_in1, i, &fe_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_1_0 = (in2->size[0] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] = in1_data[i * stride_0_0] - in2_data[i * stride_1_0];
  }
  i = in1->size[0];
  in1->size[0] = loop_ub;
  emxEnsureCapacity_real32_T(sp, in1, i, &fe_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real32_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

void b_detrend(const emlrtStack *sp, emxArray_real32_T *x)
{
  static const char_T b_fname[14] = {'L', 'A', 'P', 'A', 'C', 'K', 'E',
                                     '_', 's', 'o', 'r', 'm', 'q', 'r'};
  static const char_T fname[14] = {'L', 'A', 'P', 'A', 'C', 'K', 'E',
                                   '_', 's', 'g', 'e', 'q', 'p', '3'};
  ptrdiff_t jpvt_t[2];
  ptrdiff_t info_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  emxArray_int32_T *sOriginalType;
  emxArray_int32_T *y;
  emxArray_real32_T *A;
  emxArray_real32_T *B;
  emxArray_real32_T *W;
  emxArray_real_T *a;
  real_T *a_data;
  int32_T i;
  int32_T j;
  int32_T k;
  int32_T minmn;
  int32_T *sOriginalType_data;
  int32_T *y_data;
  real32_T b_p[2];
  real32_T tau_data[2];
  real32_T beta1;
  real32_T tol;
  real32_T *A_data;
  real32_T *B_data;
  real32_T *W_data;
  real32_T *x_data;
  char_T TRANSA1;
  char_T TRANSB1;
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
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  x_data = x->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_int32_T(sp, &y, 2, &ad_emlrtRTEI);
  st.site = &lb_emlrtRSI;
  b_st.site = &nb_emlrtRSI;
  c_st.site = &ob_emlrtRSI;
  eml_integer_colon_dispatcher(&c_st, x->size[0] - 1, y);
  y_data = y->data;
  emxInit_int32_T(sp, &sOriginalType, 1, &yd_emlrtRTEI);
  minmn = y->size[1];
  i = sOriginalType->size[0];
  sOriginalType->size[0] = y->size[1];
  emxEnsureCapacity_int32_T(sp, sOriginalType, i, &rd_emlrtRTEI);
  sOriginalType_data = sOriginalType->data;
  for (i = 0; i < minmn; i++) {
    sOriginalType_data[i] = y_data[i];
  }
  emxInit_real_T(sp, &a, 1, &ae_emlrtRTEI);
  emxInit_real32_T(sp, &W, 2, &be_emlrtRTEI);
  emxInit_real32_T(sp, &A, 2, &ce_emlrtRTEI);
  emxInit_real32_T(sp, &B, 1, &vd_emlrtRTEI);
  if (x->size[0] != 0) {
    if (x->size[0] == 1) {
      minmn = x->size[0];
      for (i = 0; i < minmn; i++) {
        x_data[i] *= 0.0F;
      }
    } else {
      int32_T jpvt[2];
      int32_T maxmn;
      boolean_T p;
      st.site = &mb_emlrtRSI;
      b_st.site = &ec_emlrtRSI;
      i = a->size[0];
      a->size[0] = y->size[1];
      emxEnsureCapacity_real_T(&b_st, a, i, &sd_emlrtRTEI);
      a_data = a->data;
      c_st.site = &fc_emlrtRSI;
      if (sOriginalType->size[0] > 2147483646) {
        d_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }
      for (i = 0; i < minmn; i++) {
        a_data[i] =
            (real_T)sOriginalType_data[i] - (real_T)sOriginalType_data[0];
      }
      for (i = 0; i < minmn; i++) {
        a_data[i] /= (real_T)sOriginalType_data[sOriginalType->size[0] - 1];
      }
      for (i = 0; i < minmn; i++) {
        real_T varargin_1;
        varargin_1 = a_data[i];
        a_data[i] = muDoubleScalarMax(varargin_1, 0.0);
      }
      i = W->size[0] * W->size[1];
      W->size[0] = y->size[1];
      W->size[1] = 2;
      emxEnsureCapacity_real32_T(&b_st, W, i, &td_emlrtRTEI);
      W_data = W->data;
      c_st.site = &gc_emlrtRSI;
      for (i = 0; i < minmn; i++) {
        c_st.site = &hc_emlrtRSI;
        d_st.site = &bb_emlrtRSI;
        e_st.site = &cb_emlrtRSI;
        W_data[i] = (real32_T)a_data[i];
      }
      c_st.site = &ic_emlrtRSI;
      for (i = 0; i < minmn; i++) {
        W_data[i + W->size[0]] = 1.0F;
      }
      c_st.site = &jc_emlrtRSI;
      d_st.site = &mc_emlrtRSI;
      i = A->size[0] * A->size[1];
      A->size[0] = y->size[1];
      A->size[1] = 2;
      emxEnsureCapacity_real32_T(&d_st, A, i, &ud_emlrtRTEI);
      A_data = A->data;
      maxmn = W->size[0] << 1;
      for (i = 0; i < maxmn; i++) {
        A_data[i] = W_data[i];
      }
      e_st.site = &oc_emlrtRSI;
      maxmn = muIntScalarMin_sint32(minmn, 2);
      if (A->size[0] == 0) {
        if (maxmn - 1 >= 0) {
          memset(&tau_data[0], 0, (uint32_T)maxmn * sizeof(real32_T));
        }
        jpvt[0] = 1;
        jpvt[1] = 2;
      } else {
        jpvt_t[0] = (ptrdiff_t)0;
        jpvt_t[1] = (ptrdiff_t)0;
        info_t =
            LAPACKE_sgeqp3(102, (ptrdiff_t)A->size[0], (ptrdiff_t)2, &A_data[0],
                           (ptrdiff_t)A->size[0], &jpvt_t[0], &tau_data[0]);
        f_st.site = &pc_emlrtRSI;
        if ((int32_T)info_t != 0) {
          p = true;
          if ((int32_T)info_t != -4) {
            if ((int32_T)info_t == -1010) {
              emlrtErrorWithMessageIdR2018a(&f_st, &w_emlrtRTEI, "MATLAB:nomem",
                                            "MATLAB:nomem", 0);
            } else {
              emlrtErrorWithMessageIdR2018a(
                  &f_st, &v_emlrtRTEI, "Coder:toolbox:LAPACKCallErrorInfo",
                  "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, &fname[0], 12,
                  (int32_T)info_t);
            }
          }
        } else {
          p = false;
        }
        if (p) {
          for (j = 0; j < 2; j++) {
            f_st.site = &qc_emlrtRSI;
            for (i = 0; i < minmn; i++) {
              A_data[j * minmn + i] = rtNaNF;
            }
          }
          minmn = maxmn - 1;
          f_st.site = &rc_emlrtRSI;
          for (k = 0; k <= minmn; k++) {
            tau_data[k] = rtNaNF;
          }
          f_st.site = &sc_emlrtRSI;
          if (maxmn + 1 <= maxmn) {
            tau_data[1] = 0.0F;
          }
          jpvt[0] = 1;
          jpvt[1] = 2;
        } else {
          jpvt[0] = (int32_T)jpvt_t[0];
          jpvt[1] = (int32_T)jpvt_t[1];
        }
      }
      k = 0;
      if (A->size[0] < 2) {
        minmn = A->size[0];
        maxmn = 2;
      } else {
        minmn = 2;
        maxmn = A->size[0];
      }
      if (minmn > 0) {
        tol = muSingleScalarMin(0.000345266977F,
                                1.1920929E-6F * (real32_T)maxmn) *
              muSingleScalarAbs(A_data[0]);
        while ((k < minmn) &&
               (!(muSingleScalarAbs(A_data[k + A->size[0] * k]) <= tol))) {
          k++;
        }
      }
      d_st.site = &nc_emlrtRSI;
      minmn = x->size[0];
      i = B->size[0];
      B->size[0] = minmn;
      emxEnsureCapacity_real32_T(&d_st, B, i, &vd_emlrtRTEI);
      B_data = B->data;
      for (i = 0; i < minmn; i++) {
        B_data[i] = x_data[i];
      }
      b_p[0] = 0.0F;
      b_p[1] = 0.0F;
      e_st.site = &tc_emlrtRSI;
      f_st.site = &wc_emlrtRSI;
      if ((A->size[0] != 0) && (B->size[0] != 0)) {
        info_t = (ptrdiff_t)B->size[0];
        info_t = LAPACKE_sormqr(102, 'L', 'T', info_t, (ptrdiff_t)1,
                                (ptrdiff_t)muIntScalarMin_sint32(A->size[0], 2),
                                &A_data[0], (ptrdiff_t)A->size[0], &tau_data[0],
                                &B_data[0], info_t);
        g_st.site = &xc_emlrtRSI;
        if ((int32_T)info_t != 0) {
          boolean_T c_p;
          p = true;
          c_p = false;
          if ((int32_T)info_t == -7) {
            c_p = true;
          } else if ((int32_T)info_t == -9) {
            c_p = true;
          } else if ((int32_T)info_t == -10) {
            c_p = true;
          }
          if (!c_p) {
            if ((int32_T)info_t == -1010) {
              emlrtErrorWithMessageIdR2018a(&g_st, &w_emlrtRTEI, "MATLAB:nomem",
                                            "MATLAB:nomem", 0);
            } else {
              emlrtErrorWithMessageIdR2018a(
                  &g_st, &v_emlrtRTEI, "Coder:toolbox:LAPACKCallErrorInfo",
                  "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, &b_fname[0],
                  12, (int32_T)info_t);
            }
          }
        } else {
          p = false;
        }
        if (p) {
          maxmn = B->size[0];
          i = B->size[0];
          B->size[0] = maxmn;
          emxEnsureCapacity_real32_T(&f_st, B, i, &wd_emlrtRTEI);
          B_data = B->data;
          for (i = 0; i < maxmn; i++) {
            B_data[i] = rtNaNF;
          }
        }
      }
      e_st.site = &uc_emlrtRSI;
      for (i = 0; i < k; i++) {
        b_p[jpvt[i] - 1] = B_data[i];
      }
      for (j = k; j >= 1; j--) {
        maxmn = jpvt[j - 1] - 1;
        b_p[maxmn] /= A_data[(j + A->size[0] * (j - 1)) - 1];
        e_st.site = &vc_emlrtRSI;
        for (i = 0; i <= j - 2; i++) {
          b_p[jpvt[0] - 1] -= b_p[maxmn] * A_data[A->size[0] * (j - 1)];
        }
      }
      if ((W->size[0] < 2) || (k < 2)) {
        c_st.site = &kc_emlrtRSI;
        warning(&c_st);
      }
      c_st.site = &lc_emlrtRSI;
      d_st.site = &yc_emlrtRSI;
      if (W->size[0] == 0) {
        B->size[0] = 0;
      } else {
        e_st.site = &ad_emlrtRSI;
        f_st.site = &bd_emlrtRSI;
        TRANSB1 = 'N';
        TRANSA1 = 'N';
        tol = 1.0F;
        beta1 = 0.0F;
        info_t = (ptrdiff_t)W->size[0];
        n_t = (ptrdiff_t)1;
        k_t = (ptrdiff_t)2;
        lda_t = (ptrdiff_t)W->size[0];
        ldb_t = (ptrdiff_t)2;
        ldc_t = (ptrdiff_t)W->size[0];
        i = B->size[0];
        B->size[0] = y->size[1];
        emxEnsureCapacity_real32_T(&f_st, B, i, &xd_emlrtRTEI);
        B_data = B->data;
        sgemm(&TRANSA1, &TRANSB1, &info_t, &n_t, &k_t, &tol, &W_data[0], &lda_t,
              &b_p[0], &ldb_t, &beta1, &B_data[0], &ldc_t);
      }
      if (x->size[0] == B->size[0]) {
        for (i = 0; i < minmn; i++) {
          x_data[i] -= B_data[i];
        }
      } else {
        c_st.site = &lc_emlrtRSI;
        minus(&c_st, x, B);
      }
    }
  }
  emxFree_real32_T(sp, &B);
  emxFree_real32_T(sp, &A);
  emxFree_real32_T(sp, &W);
  emxFree_real_T(sp, &a);
  emxFree_int32_T(sp, &y);
  emxFree_int32_T(sp, &sOriginalType);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

void detrend(const emlrtStack *sp, const emxArray_real32_T *x,
             emxArray_real32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack st;
  emxArray_int32_T *b_y;
  int32_T b;
  int32_T k;
  const real32_T *x_data;
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
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  i_st.prev = &h_st;
  i_st.tls = h_st.tls;
  x_data = x->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_int32_T(sp, &b_y, 2, &ad_emlrtRTEI);
  st.site = &lb_emlrtRSI;
  b_st.site = &nb_emlrtRSI;
  c_st.site = &ob_emlrtRSI;
  eml_integer_colon_dispatcher(&c_st, x->size[0] - 1, b_y);
  if (x->size[0] == 0) {
    y->size[0] = 0;
  } else if (x->size[0] == 1) {
    int32_T begBlock;
    begBlock = y->size[0];
    y->size[0] = 1;
    emxEnsureCapacity_real32_T(sp, y, begBlock, &xc_emlrtRTEI);
    y_data = y->data;
    y_data[0] = x_data[0] * 0.0F;
  } else {
    int32_T begBlock;
    int32_T endSeg;
    real32_T segsum;
    real32_T xs;
    st.site = &mb_emlrtRSI;
    b_st.site = &ub_emlrtRSI;
    c_st.site = &vb_emlrtRSI;
    d_st.site = &wb_emlrtRSI;
    e_st.site = &xb_emlrtRSI;
    endSeg = b_y->size[1];
    f_st.site = &yb_emlrtRSI;
    segsum = 0.0F;
    if (b_y->size[1] > 1024) {
      int32_T nblocks;
      nblocks = b_y->size[1] / 1024;
      g_st.site = &ac_emlrtRSI;
      for (b = 0; b < nblocks; b++) {
        int32_T b_endSeg;
        begBlock = b << 10;
        b_endSeg = begBlock + 1024;
        g_st.site = &bc_emlrtRSI;
        xs = 0.0F;
        h_st.site = &dc_emlrtRSI;
        for (k = begBlock + 1; k <= b_endSeg; k++) {
          xs += x_data[k - 1];
        }
        segsum += xs;
      }
      begBlock = (nblocks << 10) + 1;
    } else {
      begBlock = 1;
    }
    if (b_y->size[1] >= begBlock) {
      g_st.site = &cc_emlrtRSI;
      xs = 0.0F;
      h_st.site = &dc_emlrtRSI;
      if ((begBlock <= b_y->size[1]) && (b_y->size[1] > 2147483646)) {
        i_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }
      for (k = begBlock; k <= endSeg; k++) {
        xs += x_data[k - 1];
      }
      segsum += xs;
    }
    segsum /= (real32_T)b_y->size[1];
    begBlock = x->size[0];
    endSeg = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&b_st, y, endSeg, &yc_emlrtRTEI);
    y_data = y->data;
    for (k = 0; k < begBlock; k++) {
      y_data[k] = x_data[k] - segsum;
    }
  }
  emxFree_int32_T(sp, &b_y);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (detrend.c) */
