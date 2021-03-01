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
#include "qr.h"
#include "rt_nonfinite.h"
#include "trisolve.h"
#include "warning.h"
#include "blas.h"
#include "mwmathutil.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo qb_emlrtRSI = { 100,/* lineNo */
  "detrend",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo rb_emlrtRSI = { 174,/* lineNo */
  "detrend",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo sb_emlrtRSI = { 28, /* lineNo */
  "colon",                             /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 81, /* lineNo */
  "colon",                             /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo ac_emlrtRSI = { 184,/* lineNo */
  "chooseDetrendMethod",               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo bc_emlrtRSI = { 345,/* lineNo */
  "CCDetrend",                         /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo cc_emlrtRSI = { 63, /* lineNo */
  "applyVectorFunction",               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/applyVectorFunction.m"/* pathName */
};

static emlrtRSInfo dc_emlrtRSI = { 47, /* lineNo */
  "function_handle/parenReference",    /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/function_handle.m"/* pathName */
};

static emlrtRSInfo ec_emlrtRSI = { 344,/* lineNo */
  "@(x)subsum(x,ONE,endSeg)",          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo fc_emlrtRSI = { 244,/* lineNo */
  "subsum",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo gc_emlrtRSI = { 247,/* lineNo */
  "subsum",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo hc_emlrtRSI = { 254,/* lineNo */
  "subsum",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo ic_emlrtRSI = { 264,/* lineNo */
  "sumab",                             /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo jc_emlrtRSI = { 186,/* lineNo */
  "chooseDetrendMethod",               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo kc_emlrtRSI = { 367,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo lc_emlrtRSI = { 372,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo mc_emlrtRSI = { 377,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo nc_emlrtRSI = { 378,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo oc_emlrtRSI = { 382,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo pc_emlrtRSI = { 385,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo qc_emlrtRSI = { 387,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo rc_emlrtRSI = { 389,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo sc_emlrtRSI = { 391,/* lineNo */
  "CPPDetrend",                        /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pathName */
};

static emlrtRSInfo tc_emlrtRSI = { 14, /* lineNo */
  "max",                               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/max.m"/* pathName */
};

static emlrtRSInfo uc_emlrtRSI = { 29, /* lineNo */
  "minOrMax",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo vc_emlrtRSI = { 58, /* lineNo */
  "maximum2",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo wc_emlrtRSI = { 64, /* lineNo */
  "binaryMinOrMax",                    /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/binaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo sd_emlrtRSI = { 79, /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pathName */
};

static emlrtRSInfo ud_emlrtRSI = { 142,/* lineNo */
  "mtimes",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pathName */
};

static emlrtRSInfo vd_emlrtRSI = { 178,/* lineNo */
  "mtimes_blas",                       /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pathName */
};

static emlrtRSInfo wd_emlrtRSI = { 100,/* lineNo */
  "linsolve",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRSInfo xd_emlrtRSI = { 423,/* lineNo */
  "solveUT",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRSInfo yd_emlrtRSI = { 440,/* lineNo */
  "solveUT",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRSInfo ae_emlrtRSI = { 442,/* lineNo */
  "solveUT",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRSInfo fe_emlrtRSI = { 173,/* lineNo */
  "crudeRcondTriangular",              /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pathName */
};

static emlrtRTEInfo w_emlrtRTEI = { 88,/* lineNo */
  15,                                  /* colNo */
  "linsolve",                          /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/linsolve.m"/* pName */
};

static emlrtRTEInfo vc_emlrtRTEI = { 120,/* lineNo */
  5,                                   /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo wc_emlrtRTEI = { 124,/* lineNo */
  5,                                   /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo xc_emlrtRTEI = { 60,/* lineNo */
  20,                                  /* colNo */
  "bsxfun",                            /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/bsxfun.m"/* pName */
};

static emlrtRTEInfo yc_emlrtRTEI = { 100,/* lineNo */
  66,                                  /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo ud_emlrtRTEI = { 100,/* lineNo */
  59,                                  /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo vd_emlrtRTEI = { 365,/* lineNo */
  20,                                  /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo wd_emlrtRTEI = { 64,/* lineNo */
  67,                                  /* colNo */
  "binaryMinOrMax",                    /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/binaryMinOrMax.m"/* pName */
};

static emlrtRTEInfo xd_emlrtRTEI = { 373,/* lineNo */
  20,                                  /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo yd_emlrtRTEI = { 218,/* lineNo */
  20,                                  /* colNo */
  "mtimes",                            /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pName */
};

static emlrtRTEInfo ae_emlrtRTEI = { 100,/* lineNo */
  1,                                   /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo be_emlrtRTEI = { 365,/* lineNo */
  1,                                   /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo ce_emlrtRTEI = { 373,/* lineNo */
  1,                                   /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

static emlrtRTEInfo de_emlrtRTEI = { 1,/* lineNo */
  14,                                  /* colNo */
  "detrend",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

/* Function Definitions */
void b_detrend(const emlrtStack *sp, emxArray_real32_T *x)
{
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack j_st;
  emlrtStack st;
  emxArray_int32_T *s;
  emxArray_int32_T *y;
  emxArray_real32_T *Q;
  emxArray_real32_T *W;
  emxArray_real32_T *r;
  emxArray_real_T *a;
  emxArray_real_T *b_x;
  int32_T R_size[2];
  int32_T i;
  int32_T loop_ub;
  int32_T ns;
  int32_T nx;
  real32_T R_data[4];
  real32_T C[2];
  real32_T y_data[2];
  real32_T absAkk;
  real32_T alpha1;
  real32_T beta1;
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
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  i_st.prev = &h_st;
  i_st.tls = h_st.tls;
  j_st.prev = &i_st;
  j_st.tls = i_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_int32_T(sp, &s, 1, &ae_emlrtRTEI, true);
  emxInit_int32_T(sp, &y, 2, &yc_emlrtRTEI, true);
  st.site = &qb_emlrtRSI;
  b_st.site = &sb_emlrtRSI;
  c_st.site = &tb_emlrtRSI;
  eml_integer_colon_dispatcher(&c_st, x->size[0] - 1, y);
  i = s->size[0];
  s->size[0] = y->size[1];
  emxEnsureCapacity_int32_T(sp, s, i, &ud_emlrtRTEI);
  loop_ub = y->size[1];
  for (i = 0; i < loop_ub; i++) {
    s->data[i] = y->data[i];
  }

  emxFree_int32_T(&y);
  if (x->size[0] != 0) {
    if (x->size[0] == 1) {
      loop_ub = x->size[0];
      for (i = 0; i < loop_ub; i++) {
        x->data[i] *= 0.0F;
      }
    } else {
      emxInit_real_T(sp, &a, 1, &be_emlrtRTEI, true);
      st.site = &rb_emlrtRSI;
      b_st.site = &jc_emlrtRSI;
      ns = s->size[0] - 1;
      i = a->size[0];
      a->size[0] = s->size[0];
      emxEnsureCapacity_real_T(&b_st, a, i, &vd_emlrtRTEI);
      c_st.site = &kc_emlrtRSI;
      if ((1 <= s->size[0]) && (s->size[0] > 2147483646)) {
        d_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      for (nx = 0; nx <= ns; nx++) {
        a->data[nx] = (real_T)s->data[nx] - (real_T)s->data[0];
      }

      emxInit_real_T(&b_st, &b_x, 1, &wd_emlrtRTEI, true);
      i = b_x->size[0];
      b_x->size[0] = a->size[0];
      emxEnsureCapacity_real_T(&b_st, b_x, i, &wd_emlrtRTEI);
      loop_ub = a->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_x->data[i] = a->data[i] / (real_T)s->data[s->size[0] - 1];
      }

      c_st.site = &lc_emlrtRSI;
      d_st.site = &tc_emlrtRSI;
      e_st.site = &uc_emlrtRSI;
      f_st.site = &vc_emlrtRSI;
      g_st.site = &wc_emlrtRSI;
      i = a->size[0];
      a->size[0] = b_x->size[0];
      emxEnsureCapacity_real_T(&g_st, a, i, &id_emlrtRTEI);
      h_st.site = &xc_emlrtRSI;
      nx = b_x->size[0];
      i_st.site = &yc_emlrtRSI;
      if ((1 <= b_x->size[0]) && (b_x->size[0] > 2147483646)) {
        j_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&j_st);
      }

      for (loop_ub = 0; loop_ub < nx; loop_ub++) {
        a->data[loop_ub] = muDoubleScalarMax(b_x->data[loop_ub], 0.0);
      }

      emxFree_real_T(&b_x);
      emxInit_real32_T(&h_st, &W, 2, &ce_emlrtRTEI, true);
      i = W->size[0] * W->size[1];
      W->size[0] = s->size[0];
      W->size[1] = 2;
      emxEnsureCapacity_real32_T(&b_st, W, i, &xd_emlrtRTEI);
      c_st.site = &mc_emlrtRSI;
      if ((1 <= s->size[0]) && (s->size[0] > 2147483646)) {
        d_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      for (nx = 0; nx <= ns; nx++) {
        c_st.site = &nc_emlrtRSI;
        d_st.site = &gb_emlrtRSI;
        e_st.site = &hb_emlrtRSI;
        W->data[nx] = (real32_T)a->data[nx];
      }

      emxFree_real_T(&a);
      c_st.site = &oc_emlrtRSI;
      if ((1 <= s->size[0]) && (s->size[0] > 2147483646)) {
        d_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      for (nx = 0; nx <= ns; nx++) {
        W->data[nx + W->size[0]] = 1.0F;
      }

      emxInit_real32_T(&b_st, &Q, 2, &de_emlrtRTEI, true);
      c_st.site = &pc_emlrtRSI;
      qr(&c_st, W, Q, R_data, R_size);
      c_st.site = &qc_emlrtRSI;
      d_st.site = &td_emlrtRSI;
      if (Q->size[0] != x->size[0]) {
        if (((Q->size[0] == 1) && (Q->size[1] == 1)) || (x->size[0] == 1)) {
          emlrtErrorWithMessageIdR2018a(&d_st, &q_emlrtRTEI,
            "Coder:toolbox:mtimes_noDynamicScalarExpansion",
            "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
        } else {
          emlrtErrorWithMessageIdR2018a(&d_st, &p_emlrtRTEI, "MATLAB:innerdim",
            "MATLAB:innerdim", 0);
        }
      }

      d_st.site = &sd_emlrtRSI;
      if ((Q->size[0] == 0) || (Q->size[1] == 0) || (x->size[0] == 0)) {
        nx = Q->size[1];
        loop_ub = Q->size[1];
        if (0 <= loop_ub - 1) {
          memset(&y_data[0], 0, loop_ub * sizeof(real32_T));
        }
      } else {
        e_st.site = &ud_emlrtRSI;
        f_st.site = &vd_emlrtRSI;
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

      emxFree_real32_T(&Q);
      c_st.site = &qc_emlrtRSI;
      if (R_size[0] != nx) {
        emlrtErrorWithMessageIdR2018a(&c_st, &w_emlrtRTEI, "MATLAB:dimagree",
          "MATLAB:dimagree", 0);
      }

      d_st.site = &wd_emlrtRSI;
      nx = R_size[0];
      e_st.site = &xd_emlrtRSI;
      if (0 <= nx - 1) {
        memcpy(&C[0], &y_data[0], nx * sizeof(real32_T));
      }

      i = R_size[0] + 1;
      for (nx = i; nx < 3; nx++) {
        C[nx - 1] = 0.0F;
      }

      e_st.site = &yd_emlrtRSI;
      trisolve(R_data, R_size, C);
      e_st.site = &ae_emlrtRSI;
      if (R_size[0] == 0) {
        alpha1 = rtInfF;
      } else {
        alpha1 = muSingleScalarAbs(R_data[0]);
        beta1 = alpha1;
        nx = R_size[0];
        f_st.site = &fe_emlrtRSI;
        for (loop_ub = 2; loop_ub <= nx; loop_ub++) {
          absAkk = muSingleScalarAbs(R_data[R_size[0] + 1]);
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

      if ((2 > R_size[0]) || (alpha1 < 1.0E-10)) {
        c_st.site = &rc_emlrtRSI;
        warning(&c_st);
      }

      c_st.site = &sc_emlrtRSI;
      d_st.site = &sd_emlrtRSI;
      emxInit_real32_T(&d_st, &r, 1, &de_emlrtRTEI, true);
      if (W->size[0] == 0) {
        r->size[0] = 0;
      } else {
        e_st.site = &ud_emlrtRSI;
        f_st.site = &vd_emlrtRSI;
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
        i = r->size[0];
        r->size[0] = W->size[0];
        emxEnsureCapacity_real32_T(&f_st, r, i, &yd_emlrtRTEI);
        sgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, &W->data[0], &lda_t,
              &C[0], &ldb_t, &beta1, &r->data[0], &ldc_t);
      }

      emxFree_real32_T(&W);
      loop_ub = x->size[0];
      for (i = 0; i < loop_ub; i++) {
        x->data[i] -= r->data[i];
      }

      emxFree_real32_T(&r);
    }
  }

  emxFree_int32_T(&s);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void detrend(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T
             *y)
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
  int32_T begBlock;
  int32_T endSeg;
  int32_T k;
  int32_T nblocks;
  real32_T segsum;
  real32_T xs;
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
  emxInit_int32_T(sp, &b_y, 2, &yc_emlrtRTEI, true);
  st.site = &qb_emlrtRSI;
  b_st.site = &sb_emlrtRSI;
  c_st.site = &tb_emlrtRSI;
  eml_integer_colon_dispatcher(&c_st, x->size[0] - 1, b_y);
  if (x->size[0] == 0) {
    begBlock = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(sp, y, begBlock, &vc_emlrtRTEI);
    endSeg = x->size[0];
    for (begBlock = 0; begBlock < endSeg; begBlock++) {
      y->data[begBlock] = x->data[begBlock];
    }
  } else if (x->size[0] == 1) {
    begBlock = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(sp, y, begBlock, &wc_emlrtRTEI);
    endSeg = x->size[0];
    for (begBlock = 0; begBlock < endSeg; begBlock++) {
      y->data[begBlock] = x->data[begBlock] * 0.0F;
    }
  } else {
    st.site = &rb_emlrtRSI;
    b_st.site = &ac_emlrtRSI;
    c_st.site = &bc_emlrtRSI;
    d_st.site = &cc_emlrtRSI;
    e_st.site = &dc_emlrtRSI;
    endSeg = b_y->size[1];
    f_st.site = &ec_emlrtRSI;
    segsum = 0.0F;
    if (b_y->size[1] > 1024) {
      nblocks = b_y->size[1] / 1024;
      g_st.site = &fc_emlrtRSI;
      for (b = 0; b < nblocks; b++) {
        begBlock = (b << 10) + 1024;
        g_st.site = &gc_emlrtRSI;
        xs = 0.0F;
        h_st.site = &ic_emlrtRSI;
        if ((begBlock - 1023 <= begBlock) && (begBlock > 2147483646)) {
          i_st.site = &ab_emlrtRSI;
          check_forloop_overflow_error(&i_st);
        }

        for (k = begBlock - 1023; k <= begBlock; k++) {
          xs += x->data[k - 1];
        }

        segsum += xs;
      }

      begBlock = (nblocks << 10) + 1;
    } else {
      begBlock = 1;
    }

    if (b_y->size[1] >= begBlock) {
      g_st.site = &hc_emlrtRSI;
      xs = 0.0F;
      h_st.site = &ic_emlrtRSI;
      if ((begBlock <= b_y->size[1]) && (b_y->size[1] > 2147483646)) {
        i_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      for (k = begBlock; k <= endSeg; k++) {
        xs += x->data[k - 1];
      }

      segsum += xs;
    }

    segsum /= (real32_T)b_y->size[1];
    begBlock = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real32_T(&b_st, y, begBlock, &xc_emlrtRTEI);
    begBlock = x->size[0] - 1;
    for (k = 0; k <= begBlock; k++) {
      y->data[k] = x->data[k] - segsum;
    }
  }

  emxFree_int32_T(&b_y);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (detrend.c) */
