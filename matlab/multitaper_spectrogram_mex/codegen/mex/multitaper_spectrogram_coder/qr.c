/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * qr.c
 *
 * Code generation for function 'qr'
 *
 */

/* Include files */
#include "qr.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"
#include "lapacke.h"
#include "mwmathutil.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo ad_emlrtRSI = { 25, /* lineNo */
  "qr",                                /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/qr.m"/* pathName */
};

static emlrtRSInfo bd_emlrtRSI = { 35, /* lineNo */
  "eml_qr",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo cd_emlrtRSI = { 153,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo dd_emlrtRSI = { 161,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo ed_emlrtRSI = { 162,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo fd_emlrtRSI = { 165,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo gd_emlrtRSI = { 170,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo hd_emlrtRSI = { 174,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo id_emlrtRSI = { 176,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo jd_emlrtRSI = { 177,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo kd_emlrtRSI = { 27, /* lineNo */
  "xgeqrf",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo ld_emlrtRSI = { 91, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo md_emlrtRSI = { 94, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo nd_emlrtRSI = { 99, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo od_emlrtRSI = { 102,/* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo pd_emlrtRSI = { 189,/* lineNo */
  "unpackQR",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo qd_emlrtRSI = { 14, /* lineNo */
  "xorgqr",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

static emlrtRSInfo rd_emlrtRSI = { 60, /* lineNo */
  "ceval_xorgqr",                      /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

static emlrtRTEInfo j_emlrtRTEI = { 47,/* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

static emlrtRTEInfo k_emlrtRTEI = { 44,/* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

static emlrtRTEInfo bd_emlrtRTEI = { 35,/* lineNo */
  25,                                  /* colNo */
  "eml_qr",                            /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

static emlrtRTEInfo cd_emlrtRTEI = { 175,/* lineNo */
  20,                                  /* colNo */
  "eml_qr",                            /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

static emlrtRTEInfo dd_emlrtRTEI = { 174,/* lineNo */
  1,                                   /* colNo */
  "eml_qr",                            /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

static emlrtRTEInfo ed_emlrtRTEI = { 1,/* lineNo */
  20,                                  /* colNo */
  "qr",                                /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/matfun/qr.m"/* pName */
};

/* Function Definitions */
void qr(const emlrtStack *sp, const emxArray_real32_T *A, emxArray_real32_T *Q,
        real32_T R_data[], int32_T R_size[2])
{
  static const char_T b_fname[14] = { 'L', 'A', 'P', 'A', 'C', 'K', 'E', '_',
    's', 'o', 'r', 'g', 'q', 'r' };

  static const char_T fname[14] = { 'L', 'A', 'P', 'A', 'C', 'K', 'E', '_', 's',
    'g', 'e', 'q', 'r', 'f' };

  ptrdiff_t info_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_real32_T *b_A;
  int32_T i;
  int32_T info;
  int32_T m;
  int32_T ma;
  int32_T minmana;
  real32_T tau_data[2];
  boolean_T b_p;
  boolean_T p;
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
  emxInit_real32_T(sp, &b_A, 2, &ed_emlrtRTEI, true);
  st.site = &ad_emlrtRSI;
  b_st.site = &bd_emlrtRSI;
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = 2;
  emxEnsureCapacity_real32_T(&b_st, b_A, i, &bd_emlrtRTEI);
  ma = A->size[0] * A->size[1];
  for (i = 0; i < ma; i++) {
    b_A->data[i] = A->data[i];
  }

  c_st.site = &cd_emlrtRSI;
  m = b_A->size[0];
  d_st.site = &kd_emlrtRSI;
  ma = b_A->size[0];
  minmana = muIntScalarMin_sint32(ma, 2);
  if (b_A->size[0] == 0) {
    if (0 <= minmana - 1) {
      memset(&tau_data[0], 0, minmana * sizeof(real32_T));
    }
  } else {
    info_t = LAPACKE_sgeqrf(102, (ptrdiff_t)b_A->size[0], (ptrdiff_t)2,
      &b_A->data[0], (ptrdiff_t)b_A->size[0], &tau_data[0]);
    info = (int32_T)info_t;
    e_st.site = &ld_emlrtRSI;
    if (info != 0) {
      p = true;
      if (info != -4) {
        if (info == -1010) {
          emlrtErrorWithMessageIdR2018a(&e_st, &k_emlrtRTEI, "MATLAB:nomem",
            "MATLAB:nomem", 0);
        } else {
          emlrtErrorWithMessageIdR2018a(&e_st, &j_emlrtRTEI,
            "Coder:toolbox:LAPACKCallErrorInfo",
            "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, fname, 12, info);
        }
      }
    } else {
      p = false;
    }

    if (p) {
      e_st.site = &md_emlrtRSI;
      if ((1 <= m) && (m > 2147483646)) {
        f_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      for (info = 0; info < m; info++) {
        b_A->data[info] = rtNaNF;
      }

      e_st.site = &md_emlrtRSI;
      for (info = 0; info < m; info++) {
        b_A->data[ma + info] = rtNaNF;
      }

      ma = muIntScalarMin_sint32(m, 2) - 1;
      e_st.site = &nd_emlrtRSI;
      for (info = 0; info <= ma; info++) {
        tau_data[info] = rtNaNF;
      }

      e_st.site = &od_emlrtRSI;
      if (ma + 2 <= minmana) {
        tau_data[1] = 0.0F;
      }
    }
  }

  m = b_A->size[0];
  minmana = muIntScalarMin_sint32(m, 2);
  R_size[0] = minmana;
  R_size[1] = 2;
  c_st.site = &dd_emlrtRSI;
  for (ma = 0; ma < minmana; ma++) {
    c_st.site = &ed_emlrtRSI;
    for (info = 0; info <= ma; info++) {
      R_data[info + R_size[0] * ma] = b_A->data[info + b_A->size[0] * ma];
    }

    c_st.site = &fd_emlrtRSI;
    if (ma + 2 <= minmana) {
      R_data[R_size[0] * ma + 1] = 0.0F;
    }
  }

  i = b_A->size[0] + 1;
  for (ma = i; ma < 3; ma++) {
    c_st.site = &gd_emlrtRSI;
    for (info = 0; info < minmana; info++) {
      R_data[info + R_size[0] * (ma - 1)] = b_A->data[info + b_A->size[0] * (ma
        - 1)];
    }
  }

  c_st.site = &hd_emlrtRSI;
  d_st.site = &pd_emlrtRSI;
  e_st.site = &qd_emlrtRSI;
  if (b_A->size[0] != 0) {
    info_t = LAPACKE_sorgqr(102, (ptrdiff_t)b_A->size[0], (ptrdiff_t)minmana,
      (ptrdiff_t)minmana, &b_A->data[0], (ptrdiff_t)b_A->size[0], &tau_data[0]);
    info = (int32_T)info_t;
    f_st.site = &rd_emlrtRSI;
    if (info != 0) {
      p = true;
      b_p = false;
      if (info == -7) {
        b_p = true;
      } else {
        if (info == -5) {
          b_p = true;
        }
      }

      if (!b_p) {
        if (info == -1010) {
          emlrtErrorWithMessageIdR2018a(&f_st, &k_emlrtRTEI, "MATLAB:nomem",
            "MATLAB:nomem", 0);
        } else {
          emlrtErrorWithMessageIdR2018a(&f_st, &j_emlrtRTEI,
            "Coder:toolbox:LAPACKCallErrorInfo",
            "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, b_fname, 12, info);
        }
      }
    } else {
      p = false;
    }

    if (p) {
      ma = b_A->size[0];
      info = b_A->size[1];
      i = b_A->size[0] * b_A->size[1];
      b_A->size[0] = ma;
      b_A->size[1] = 2;
      emxEnsureCapacity_real32_T(&e_st, b_A, i, &dd_emlrtRTEI);
      ma *= info;
      for (i = 0; i < ma; i++) {
        b_A->data[i] = rtNaNF;
      }
    }
  }

  i = Q->size[0] * Q->size[1];
  Q->size[0] = m;
  Q->size[1] = minmana;
  emxEnsureCapacity_real32_T(&b_st, Q, i, &cd_emlrtRTEI);
  c_st.site = &id_emlrtRSI;
  for (ma = 0; ma < minmana; ma++) {
    c_st.site = &jd_emlrtRSI;
    if ((1 <= m) && (m > 2147483646)) {
      d_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&d_st);
    }

    for (info = 0; info < m; info++) {
      Q->data[info + Q->size[0] * ma] = b_A->data[info + b_A->size[0] * ma];
    }
  }

  emxFree_real32_T(&b_A);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (qr.c) */
