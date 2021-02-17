/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * detrend.cpp
 *
 * Code generation for function 'detrend'
 *
 */

/* Include files */
#include "detrend.h"
#include "blas.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "mwmathutil.h"
#include "qr.h"
#include "rt_nonfinite.h"
#include "trisolve.h"
#include "warning.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo pb_emlrtRSI = { 30, /* lineNo */
  "detrend",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo qb_emlrtRSI = { 88, /* lineNo */
  "detrend",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo rb_emlrtRSI = { 130,/* lineNo */
  "detrend",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo sb_emlrtRSI = { 28, /* lineNo */
  "colon",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 81, /* lineNo */
  "colon",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo ub_emlrtRSI = { 126,/* lineNo */
  "eml_integer_colon_dispatcher",      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo vb_emlrtRSI = { 149,/* lineNo */
  "eml_signed_integer_colon",          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo wb_emlrtRSI = { 154,/* lineNo */
  "eml_signed_integer_colon",          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo xb_emlrtRSI = { 239,/* lineNo */
  "integer_colon_length",              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo yb_emlrtRSI = { 273,/* lineNo */
  "integer_colon_length_nonnegd",      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo ac_emlrtRSI = { 290,/* lineNo */
  "CCDetrend",                         /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo bc_emlrtRSI = { 63, /* lineNo */
  "applyVectorFunction",               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/applyVectorFunction.m"/* pathName */
};

static emlrtRSInfo cc_emlrtRSI = { 39, /* lineNo */
  "function_handle/parenReference",    /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/function_handle.m"/* pathName */
};

static emlrtRSInfo dc_emlrtRSI = { 289,/* lineNo */
  "@(x)subsum(x,ONE,endSeg)",          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo ec_emlrtRSI = { 190,/* lineNo */
  "subsum",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo fc_emlrtRSI = { 193,/* lineNo */
  "subsum",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo gc_emlrtRSI = { 200,/* lineNo */
  "subsum",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo hc_emlrtRSI = { 210,/* lineNo */
  "sumab",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo ic_emlrtRSI = { 132,/* lineNo */
  "detrend",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo jc_emlrtRSI = { 312,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo kc_emlrtRSI = { 317,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo lc_emlrtRSI = { 322,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo mc_emlrtRSI = { 323,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo nc_emlrtRSI = { 327,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo oc_emlrtRSI = { 330,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo pc_emlrtRSI = { 332,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo qc_emlrtRSI = { 334,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo rc_emlrtRSI = { 336,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo sc_emlrtRSI = { 14, /* lineNo */
  "max",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/max.m"/* pathName */
};

static emlrtRSInfo tc_emlrtRSI = { 20, /* lineNo */
  "minOrMax",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo uc_emlrtRSI = { 38, /* lineNo */
  "unaryOrBinaryDispatch",             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo vc_emlrtRSI = { 62, /* lineNo */
  "binaryMinOrMax",                    /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/binaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo ge_emlrtRSI = { 232,/* lineNo */
  "mtimes",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pathName */
};

static emlrtRSInfo he_emlrtRSI = { 228,/* lineNo */
  "mtimes",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pathName */
};

static emlrtRSInfo ie_emlrtRSI = { 100,/* lineNo */
  "linsolve",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRSInfo je_emlrtRSI = { 423,/* lineNo */
  "solveUT",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRSInfo ke_emlrtRSI = { 440,/* lineNo */
  "solveUT",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRSInfo le_emlrtRSI = { 442,/* lineNo */
  "solveUT",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRSInfo qe_emlrtRSI = { 173,/* lineNo */
  "crudeRcondTriangular",              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRTEInfo i_emlrtRTEI = { 88,/* lineNo */
  15,                                  /* colNo */
  "linsolve",                          /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pName */
};

static emlrtRTEInfo ad_emlrtRTEI = { 30,/* lineNo */
  18,                                  /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo bd_emlrtRTEI = { 126,/* lineNo */
  9,                                   /* colNo */
  "colon",                             /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pName */
};

static emlrtRTEInfo cd_emlrtRTEI = { 30,/* lineNo */
  5,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo dd_emlrtRTEI = { 88,/* lineNo */
  59,                                  /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo ed_emlrtRTEI = { 31,/* lineNo */
  5,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo fd_emlrtRTEI = { 292,/* lineNo */
  5,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo gd_emlrtRTEI = { 88,/* lineNo */
  66,                                  /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo hd_emlrtRTEI = { 88,/* lineNo */
  1,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo id_emlrtRTEI = { 296,/* lineNo */
  14,                                  /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo jd_emlrtRTEI = { 62,/* lineNo */
  67,                                  /* colNo */
  "binaryMinOrMax",                    /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/binaryMinOrMax.m"/* pName */
};

static emlrtRTEInfo kd_emlrtRTEI = { 62,/* lineNo */
  10,                                  /* colNo */
  "binaryMinOrMax",                    /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/binaryMinOrMax.m"/* pName */
};

static emlrtRTEInfo ld_emlrtRTEI = { 232,/* lineNo */
  13,                                  /* colNo */
  "mtimes",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pName */
};

static emlrtRTEInfo md_emlrtRTEI = { 336,/* lineNo */
  1,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo nd_emlrtRTEI = { 310,/* lineNo */
  1,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo od_emlrtRTEI = { 318,/* lineNo */
  1,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo pd_emlrtRTEI = { 330,/* lineNo */
  2,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo qd_emlrtRTEI = { 330,/* lineNo */
  4,                                   /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

/* Function Declarations */
static void CPPDetrend(const emlrtStack *sp, const emxArray_real32_T *x, real_T
  bp, const emxArray_real_T *s, emxArray_real32_T *y);

/* Function Definitions */
static void CPPDetrend(const emlrtStack *sp, const emxArray_real32_T *x, real_T
  bp, const emxArray_real_T *s, emxArray_real32_T *y)
{
  emxArray_real_T *a;
  int32_T ns;
  int32_T i;
  int32_T nx;
  emxArray_real_T *b_x;
  int32_T loop_ub;
  int32_T k;
  emxArray_real32_T *W;
  emxArray_real32_T *Q;
  emxArray_real32_T *R;
  real32_T y_data[2];
  char_T TRANSB1;
  char_T TRANSA1;
  real32_T alpha1;
  real32_T beta1;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  real32_T C[2];
  real32_T absAkk;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
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
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &a, 1, &nd_emlrtRTEI, true);
  ns = s->size[0] - 1;
  i = a->size[0];
  a->size[0] = s->size[0];
  emxEnsureCapacity_real_T(sp, a, i, &id_emlrtRTEI);
  st.site = &jc_emlrtRSI;
  if ((1 <= s->size[0]) && (s->size[0] > 2147483646)) {
    b_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }

  for (nx = 0; nx <= ns; nx++) {
    a->data[nx] = s->data[nx] - bp;
  }

  emxInit_real_T(sp, &b_x, 1, &jd_emlrtRTEI, true);
  i = b_x->size[0];
  b_x->size[0] = a->size[0];
  emxEnsureCapacity_real_T(sp, b_x, i, &jd_emlrtRTEI);
  loop_ub = a->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_x->data[i] = a->data[i] / s->data[s->size[0] - 1];
  }

  st.site = &kc_emlrtRSI;
  b_st.site = &sc_emlrtRSI;
  c_st.site = &tc_emlrtRSI;
  d_st.site = &uc_emlrtRSI;
  e_st.site = &vc_emlrtRSI;
  i = a->size[0];
  a->size[0] = b_x->size[0];
  emxEnsureCapacity_real_T(&e_st, a, i, &kd_emlrtRTEI);
  f_st.site = &wc_emlrtRSI;
  nx = b_x->size[0];
  g_st.site = &xc_emlrtRSI;
  if ((1 <= b_x->size[0]) && (b_x->size[0] > 2147483646)) {
    h_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&h_st);
  }

  for (k = 0; k < nx; k++) {
    a->data[k] = muDoubleScalarMax(b_x->data[k], 0.0);
  }

  emxFree_real_T(&b_x);
  emxInit_real32_T(&f_st, &W, 2, &od_emlrtRTEI, true);
  i = W->size[0] * W->size[1];
  W->size[0] = s->size[0];
  W->size[1] = 2;
  emxEnsureCapacity_real32_T(sp, W, i, &id_emlrtRTEI);
  st.site = &lc_emlrtRSI;
  if ((1 <= s->size[0]) && (s->size[0] > 2147483646)) {
    b_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }

  for (nx = 0; nx <= ns; nx++) {
    st.site = &mc_emlrtRSI;
    b_st.site = &gb_emlrtRSI;
    c_st.site = &hb_emlrtRSI;
    W->data[nx] = static_cast<real32_T>(a->data[nx]);
  }

  emxFree_real_T(&a);
  st.site = &nc_emlrtRSI;
  if ((1 <= s->size[0]) && (s->size[0] > 2147483646)) {
    b_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }

  for (nx = 0; nx <= ns; nx++) {
    W->data[nx + W->size[0]] = 1.0F;
  }

  emxInit_real32_T(sp, &Q, 2, &pd_emlrtRTEI, true);
  emxInit_real32_T(sp, &R, 2, &qd_emlrtRTEI, true);
  st.site = &oc_emlrtRSI;
  qr(&st, W, Q, R);
  st.site = &pc_emlrtRSI;
  b_st.site = &fe_emlrtRSI;
  if (Q->size[0] != x->size[0]) {
    if (((Q->size[0] == 1) && (Q->size[1] == 1)) || (x->size[0] == 1)) {
      emlrtErrorWithMessageIdR2018a(&b_st, &j_emlrtRTEI,
        "Coder:toolbox:mtimes_noDynamicScalarExpansion",
        "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
    } else {
      emlrtErrorWithMessageIdR2018a(&b_st, &k_emlrtRTEI, "MATLAB:innerdim",
        "MATLAB:innerdim", 0);
    }
  }

  if ((Q->size[0] == 1) || (x->size[0] == 1)) {
    nx = Q->size[1];
    loop_ub = Q->size[1];
    for (i = 0; i < loop_ub; i++) {
      y_data[i] = 0.0F;
      k = Q->size[0];
      for (ns = 0; ns < k; ns++) {
        y_data[i] += Q->data[ns + Q->size[0] * i] * x->data[ns];
      }
    }
  } else {
    b_st.site = &ee_emlrtRSI;
    if ((Q->size[0] == 0) || (Q->size[1] == 0) || (x->size[0] == 0)) {
      nx = Q->size[1];
      loop_ub = Q->size[1];
      if (0 <= loop_ub - 1) {
        memset(&y_data[0], 0, loop_ub * sizeof(real32_T));
      }
    } else {
      c_st.site = &he_emlrtRSI;
      c_st.site = &ge_emlrtRSI;
      TRANSB1 = 'N';
      TRANSA1 = 'T';
      alpha1 = 1.0F;
      beta1 = 0.0F;
      m_t = (ptrdiff_t)Q->size[1];
      n_t = (ptrdiff_t)1;
      k_t = (ptrdiff_t)Q->size[0];
      lda_t = (ptrdiff_t)Q->size[0];
      ldb_t = (ptrdiff_t)x->size[0];
      ldc_t = (ptrdiff_t)Q->size[1];
      nx = Q->size[1];
      sgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, &Q->data[0], &lda_t,
            &x->data[0], &ldb_t, &beta1, &y_data[0], &ldc_t);
    }
  }

  emxFree_real32_T(&Q);
  st.site = &pc_emlrtRSI;
  if (R->size[0] != nx) {
    emlrtErrorWithMessageIdR2018a(&st, &i_emlrtRTEI, "MATLAB:dimagree",
      "MATLAB:dimagree", 0);
  }

  b_st.site = &ie_emlrtRSI;
  nx = R->size[0];
  nx = muIntScalarMin_sint32(nx, 2) - 1;
  c_st.site = &je_emlrtRSI;
  if (0 <= nx) {
    memcpy(&C[0], &y_data[0], (nx + 1) * sizeof(real32_T));
  }

  i = nx + 2;
  for (nx = i; nx < 3; nx++) {
    C[nx - 1] = 0.0F;
  }

  c_st.site = &ke_emlrtRSI;
  trisolve(R, C);
  c_st.site = &le_emlrtRSI;
  if (R->size[0] == 0) {
    alpha1 = rtInfF;
  } else {
    alpha1 = muSingleScalarAbs(R->data[0]);
    beta1 = alpha1;
    nx = static_cast<int32_T>(muDoubleScalarMin(static_cast<real_T>(R->size[0]),
      2.0));
    d_st.site = &qe_emlrtRSI;
    for (k = 2; k <= nx; k++) {
      absAkk = muSingleScalarAbs(R->data[R->size[0] + 1]);
      if ((absAkk > alpha1) || muSingleScalarIsNaN(absAkk)) {
        alpha1 = absAkk;
      } else {
        if (absAkk < beta1) {
          beta1 = absAkk;
        }
      }
    }

    alpha1 = beta1 / alpha1;
  }

  if ((2 > R->size[0]) || (alpha1 < 1.0E-10)) {
    st.site = &qc_emlrtRSI;
    warning(&st);
  }

  emxFree_real32_T(&R);
  st.site = &rc_emlrtRSI;
  b_st.site = &ee_emlrtRSI;
  if (W->size[0] != 0) {
    c_st.site = &he_emlrtRSI;
    c_st.site = &ge_emlrtRSI;
    TRANSB1 = 'N';
    TRANSA1 = 'N';
    alpha1 = 1.0F;
    beta1 = 0.0F;
    m_t = (ptrdiff_t)W->size[0];
    n_t = (ptrdiff_t)1;
    k_t = (ptrdiff_t)2;
    lda_t = (ptrdiff_t)W->size[0];
    ldb_t = (ptrdiff_t)2;
    ldc_t = (ptrdiff_t)W->size[0];
    i = y->size[0];
    y->size[0] = W->size[0];
    emxEnsureCapacity_real32_T(&c_st, y, i, &ld_emlrtRTEI);
    sgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, &W->data[0], &lda_t,
          &C[0], &ldb_t, &beta1, &y->data[0], &ldc_t);
  }

  emxFree_real32_T(&W);
  i = y->size[0];
  y->size[0] = x->size[0];
  emxEnsureCapacity_real32_T(sp, y, i, &md_emlrtRTEI);
  loop_ub = x->size[0];
  for (i = 0; i < loop_ub; i++) {
    y->data[i] = x->data[i] - y->data[i];
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void b_detrend(const emlrtStack *sp, const emxArray_real32_T *x,
               emxArray_real32_T *y)
{
  emxArray_real32_T *b_x;
  int32_T yk;
  int32_T k;
  emxArray_int32_T *b_y;
  int32_T n;
  emxArray_real_T *s;
  emxArray_real32_T *b_y1;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
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
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real32_T(sp, &b_x, 1, &ad_emlrtRTEI, true);
  st.site = &pb_emlrtRSI;
  yk = x->size[1];
  k = b_x->size[0];
  b_x->size[0] = x->size[1];
  emxEnsureCapacity_real32_T(&st, b_x, k, &ad_emlrtRTEI);
  for (k = 0; k < yk; k++) {
    b_x->data[k] = x->data[k];
  }

  emxInit_int32_T(&st, &b_y, 2, &gd_emlrtRTEI, true);
  b_st.site = &qb_emlrtRSI;
  c_st.site = &sb_emlrtRSI;
  d_st.site = &tb_emlrtRSI;
  e_st.site = &ub_emlrtRSI;
  f_st.site = &vb_emlrtRSI;
  if (b_x->size[0] - 1 < 0) {
    n = 0;
  } else {
    n = b_x->size[0];
  }

  k = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n;
  emxEnsureCapacity_int32_T(&e_st, b_y, k, &bd_emlrtRTEI);
  if (n > 0) {
    b_y->data[0] = 0;
    yk = 0;
    f_st.site = &wb_emlrtRSI;
    if ((2 <= n) && (n > 2147483646)) {
      g_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&g_st);
    }

    for (k = 2; k <= n; k++) {
      yk++;
      b_y->data[k - 1] = yk;
    }
  }

  emxInit_real_T(&e_st, &s, 1, &hd_emlrtRTEI, true);
  k = s->size[0];
  s->size[0] = b_y->size[1];
  emxEnsureCapacity_real_T(&st, s, k, &dd_emlrtRTEI);
  yk = b_y->size[1];
  for (k = 0; k < yk; k++) {
    s->data[k] = b_y->data[k];
  }

  emxFree_int32_T(&b_y);
  emxInit_real32_T(&st, &b_y1, 1, &cd_emlrtRTEI, true);
  if (b_x->size[0] == 0) {
    k = b_y1->size[0];
    b_y1->size[0] = b_x->size[0];
    emxEnsureCapacity_real32_T(&st, b_y1, k, &cd_emlrtRTEI);
    yk = b_x->size[0];
    for (k = 0; k < yk; k++) {
      b_y1->data[k] = b_x->data[k];
    }
  } else if (b_x->size[0] == 1) {
    k = b_y1->size[0];
    b_y1->size[0] = b_x->size[0];
    emxEnsureCapacity_real32_T(&st, b_y1, k, &cd_emlrtRTEI);
    yk = b_x->size[0];
    for (k = 0; k < yk; k++) {
      b_y1->data[k] = b_x->data[k] * 0.0F;
    }
  } else {
    b_st.site = &ic_emlrtRSI;
    CPPDetrend(&b_st, b_x, s->data[0], s, b_y1);
  }

  emxFree_real_T(&s);
  emxFree_real32_T(&b_x);
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = b_y1->size[0];
  emxEnsureCapacity_real32_T(sp, y, k, &ed_emlrtRTEI);
  yk = b_y1->size[0];
  for (k = 0; k < yk; k++) {
    y->data[k] = b_y1->data[k];
  }

  emxFree_real32_T(&b_y1);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void detrend(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T
             *y)
{
  emxArray_real32_T *b_x;
  int32_T yk;
  int32_T nblocks;
  emxArray_int32_T *b_y;
  int32_T n;
  emxArray_real32_T *b_y1;
  int32_T k;
  emxArray_real_T *c_y;
  int32_T endSeg;
  real32_T segsum;
  int32_T b;
  real32_T xs;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
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
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real32_T(sp, &b_x, 1, &ad_emlrtRTEI, true);
  st.site = &pb_emlrtRSI;
  yk = x->size[1];
  nblocks = b_x->size[0];
  b_x->size[0] = x->size[1];
  emxEnsureCapacity_real32_T(&st, b_x, nblocks, &ad_emlrtRTEI);
  for (nblocks = 0; nblocks < yk; nblocks++) {
    b_x->data[nblocks] = x->data[nblocks];
  }

  emxInit_int32_T(&st, &b_y, 2, &gd_emlrtRTEI, true);
  b_st.site = &qb_emlrtRSI;
  c_st.site = &sb_emlrtRSI;
  d_st.site = &tb_emlrtRSI;
  e_st.site = &ub_emlrtRSI;
  f_st.site = &vb_emlrtRSI;
  g_st.site = &xb_emlrtRSI;
  if (b_x->size[0] - 1 < 0) {
    n = 0;
  } else {
    h_st.site = &yb_emlrtRSI;
    n = b_x->size[0];
  }

  nblocks = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n;
  emxEnsureCapacity_int32_T(&e_st, b_y, nblocks, &bd_emlrtRTEI);
  if (n > 0) {
    b_y->data[0] = 0;
    yk = 0;
    f_st.site = &wb_emlrtRSI;
    if ((2 <= n) && (n > 2147483646)) {
      g_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&g_st);
    }

    for (k = 2; k <= n; k++) {
      yk++;
      b_y->data[k - 1] = yk;
    }
  }

  emxInit_real32_T(&st, &b_y1, 1, &cd_emlrtRTEI, true);
  if (b_x->size[0] == 0) {
    nblocks = b_y1->size[0];
    b_y1->size[0] = b_x->size[0];
    emxEnsureCapacity_real32_T(&st, b_y1, nblocks, &cd_emlrtRTEI);
    yk = b_x->size[0];
    for (nblocks = 0; nblocks < yk; nblocks++) {
      b_y1->data[nblocks] = b_x->data[nblocks];
    }
  } else if (b_x->size[0] == 1) {
    nblocks = b_y1->size[0];
    b_y1->size[0] = b_x->size[0];
    emxEnsureCapacity_real32_T(&st, b_y1, nblocks, &cd_emlrtRTEI);
    yk = b_x->size[0];
    for (nblocks = 0; nblocks < yk; nblocks++) {
      b_y1->data[nblocks] = b_x->data[nblocks] * 0.0F;
    }
  } else {
    emxInit_real_T(&st, &c_y, 1, &dd_emlrtRTEI, true);
    b_st.site = &rb_emlrtRSI;
    nblocks = c_y->size[0];
    c_y->size[0] = b_y->size[1];
    emxEnsureCapacity_real_T(&b_st, c_y, nblocks, &dd_emlrtRTEI);
    yk = b_y->size[1];
    for (nblocks = 0; nblocks < yk; nblocks++) {
      c_y->data[nblocks] = b_y->data[nblocks];
    }

    endSeg = c_y->size[0];
    c_st.site = &ac_emlrtRSI;
    d_st.site = &bc_emlrtRSI;
    e_st.site = &cc_emlrtRSI;
    nblocks = c_y->size[0];
    c_y->size[0] = b_y->size[1];
    emxEnsureCapacity_real_T(&e_st, c_y, nblocks, &dd_emlrtRTEI);
    yk = b_y->size[1];
    for (nblocks = 0; nblocks < yk; nblocks++) {
      c_y->data[nblocks] = b_y->data[nblocks];
    }

    yk = c_y->size[0];
    f_st.site = &dc_emlrtRSI;
    segsum = 0.0F;
    emxFree_real_T(&c_y);
    if (endSeg > 1024) {
      nblocks = endSeg / 1024;
      g_st.site = &ec_emlrtRSI;
      for (b = 0; b < nblocks; b++) {
        n = (b << 10) + 1024;
        g_st.site = &fc_emlrtRSI;
        xs = 0.0F;
        h_st.site = &hc_emlrtRSI;
        if ((n - 1023 <= n) && (n > 2147483646)) {
          i_st.site = &ab_emlrtRSI;
          check_forloop_overflow_error(&i_st);
        }

        for (k = n - 1023; k <= n; k++) {
          xs += b_x->data[k - 1];
        }

        segsum += xs;
      }

      n = (nblocks << 10) + 1;
    } else {
      n = 1;
    }

    if (endSeg >= n) {
      g_st.site = &gc_emlrtRSI;
      xs = 0.0F;
      h_st.site = &hc_emlrtRSI;
      if ((n <= endSeg) && (endSeg > 2147483646)) {
        i_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      for (k = n; k <= yk; k++) {
        xs += b_x->data[k - 1];
      }

      segsum += xs;
    }

    segsum /= static_cast<real32_T>(endSeg);
    nblocks = b_y1->size[0];
    b_y1->size[0] = b_x->size[0];
    emxEnsureCapacity_real32_T(&b_st, b_y1, nblocks, &fd_emlrtRTEI);
    nblocks = b_x->size[0] - 1;
    for (k = 0; k <= nblocks; k++) {
      b_y1->data[k] = b_x->data[k] - segsum;
    }
  }

  emxFree_int32_T(&b_y);
  emxFree_real32_T(&b_x);
  nblocks = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = b_y1->size[0];
  emxEnsureCapacity_real32_T(sp, y, nblocks, &ed_emlrtRTEI);
  yk = b_y1->size[0];
  for (nblocks = 0; nblocks < yk; nblocks++) {
    y->data[nblocks] = b_y1->data[nblocks];
  }

  emxFree_real32_T(&b_y1);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (detrend.cpp) */
