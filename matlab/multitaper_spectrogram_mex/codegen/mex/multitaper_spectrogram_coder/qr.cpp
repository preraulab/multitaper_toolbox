/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * qr.cpp
 *
 * Code generation for function 'qr'
 *
 */

/* Include files */
#include "qr.h"
#include "eml_int_forloop_overflow_check.h"
#include "lapacke.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "rt_nonfinite.h"
#include "xgeqrf.h"

/* Variable Definitions */
static emlrtRSInfo yc_emlrtRSI = { 25, /* lineNo */
  "qr",                                /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/qr.m"/* pathName */
};

static emlrtRSInfo ad_emlrtRSI = { 35, /* lineNo */
  "eml_qr",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo bd_emlrtRSI = { 153,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo cd_emlrtRSI = { 161,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo dd_emlrtRSI = { 168,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo ed_emlrtRSI = { 173,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo fd_emlrtRSI = { 174,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo gd_emlrtRSI = { 177,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo hd_emlrtRSI = { 182,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo id_emlrtRSI = { 186,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo jd_emlrtRSI = { 188,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo kd_emlrtRSI = { 189,/* lineNo */
  "qr_econ",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo xd_emlrtRSI = { 202,/* lineNo */
  "unpackQR",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pathName */
};

static emlrtRSInfo yd_emlrtRSI = { 14, /* lineNo */
  "xorgqr",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

static emlrtRSInfo de_emlrtRSI = { 59, /* lineNo */
  "ceval_xorgqr",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

static emlrtRTEInfo rd_emlrtRTEI = { 153,/* lineNo */
  6,                                   /* colNo */
  "eml_qr",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

static emlrtRTEInfo sd_emlrtRTEI = { 35,/* lineNo */
  9,                                   /* colNo */
  "eml_qr",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

static emlrtRTEInfo td_emlrtRTEI = { 186,/* lineNo */
  5,                                   /* colNo */
  "eml_qr",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

static emlrtRTEInfo ud_emlrtRTEI = { 168,/* lineNo */
  5,                                   /* colNo */
  "eml_qr",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

static emlrtRTEInfo vd_emlrtRTEI = { 169,/* lineNo */
  5,                                   /* colNo */
  "eml_qr",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

static emlrtRTEInfo wd_emlrtRTEI = { 35,/* lineNo */
  25,                                  /* colNo */
  "eml_qr",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/private/eml_qr.m"/* pName */
};

/* Function Definitions */
void qr(const emlrtStack *sp, const emxArray_real32_T *A, emxArray_real32_T *Q,
        emxArray_real32_T *R)
{
  emxArray_real32_T *b_A;
  int32_T i;
  int32_T loop_ub;
  real32_T tau_data[2];
  int32_T tau_size[1];
  int32_T m;
  emxArray_real32_T *c_A;
  int32_T b_i;
  ptrdiff_t info_t;
  boolean_T p;
  boolean_T b_p;
  static const char_T fname[14] = { 'L', 'A', 'P', 'A', 'C', 'K', 'E', '_', 's',
    'o', 'r', 'g', 'q', 'r' };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
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
  emxInit_real32_T(sp, &b_A, 2, &rd_emlrtRTEI, true);
  st.site = &yc_emlrtRSI;
  b_st.site = &ad_emlrtRSI;
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = 2;
  emxEnsureCapacity_real32_T(&b_st, b_A, i, &rd_emlrtRTEI);
  loop_ub = A->size[0] * A->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_A->data[i] = A->data[i];
  }

  c_st.site = &bd_emlrtRSI;
  xgeqrf(&c_st, b_A, tau_data, tau_size);
  m = b_A->size[0];
  emxInit_real32_T(&b_st, &c_A, 2, &wd_emlrtRTEI, true);
  if (b_A->size[0] > 2) {
    i = R->size[0] * R->size[1];
    R->size[0] = 2;
    R->size[1] = 2;
    emxEnsureCapacity_real32_T(&b_st, R, i, &sd_emlrtRTEI);
    c_st.site = &cd_emlrtRSI;
    for (b_i = 0; b_i < 1; b_i++) {
      R->data[0] = b_A->data[0];
    }

    R->data[1] = 0.0F;
    c_st.site = &cd_emlrtRSI;
    for (b_i = 0; b_i < 2; b_i++) {
      R->data[b_i + R->size[0]] = b_A->data[b_i + b_A->size[0]];
    }

    c_st.site = &dd_emlrtRSI;
    i = c_A->size[0] * c_A->size[1];
    c_A->size[0] = b_A->size[0];
    c_A->size[1] = 2;
    emxEnsureCapacity_real32_T(&c_st, c_A, i, &ud_emlrtRTEI);
    loop_ub = b_A->size[0] * b_A->size[1];
    for (i = 0; i < loop_ub; i++) {
      c_A->data[i] = b_A->data[i];
    }

    d_st.site = &xd_emlrtRSI;
    e_st.site = &yd_emlrtRSI;
    info_t = LAPACKE_sorgqr(102, (ptrdiff_t)b_A->size[0], (ptrdiff_t)2,
      (ptrdiff_t)2, &c_A->data[0], (ptrdiff_t)b_A->size[0], &tau_data[0]);
    loop_ub = (int32_T)info_t;
    f_st.site = &de_emlrtRSI;
    if (loop_ub != 0) {
      p = true;
      b_p = false;
      if (loop_ub == -7) {
        b_p = true;
      } else {
        if (loop_ub == -5) {
          b_p = true;
        }
      }

      if (!b_p) {
        if (loop_ub == -1010) {
          emlrtErrorWithMessageIdR2018a(&f_st, &m_emlrtRTEI, "MATLAB:nomem",
            "MATLAB:nomem", 0);
        } else {
          emlrtErrorWithMessageIdR2018a(&f_st, &l_emlrtRTEI,
            "Coder:toolbox:LAPACKCallErrorInfo",
            "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, fname, 12, loop_ub);
        }
      }
    } else {
      p = false;
    }

    if (p) {
      i = c_A->size[0] * c_A->size[1];
      c_A->size[1] = 2;
      emxEnsureCapacity_real32_T(&e_st, c_A, i, &ud_emlrtRTEI);
      loop_ub = c_A->size[0];
      for (i = 0; i < loop_ub; i++) {
        c_A->data[i] = rtNaNF;
      }

      loop_ub = c_A->size[0];
      for (i = 0; i < loop_ub; i++) {
        c_A->data[i + c_A->size[0]] = rtNaNF;
      }
    }

    i = Q->size[0] * Q->size[1];
    Q->size[0] = c_A->size[0];
    Q->size[1] = 2;
    emxEnsureCapacity_real32_T(&b_st, Q, i, &vd_emlrtRTEI);
    loop_ub = c_A->size[0] * c_A->size[1];
    for (i = 0; i < loop_ub; i++) {
      Q->data[i] = c_A->data[i];
    }
  } else {
    i = R->size[0] * R->size[1];
    R->size[0] = b_A->size[0];
    R->size[1] = b_A->size[1];
    emxEnsureCapacity_real32_T(&b_st, R, i, &sd_emlrtRTEI);
    c_st.site = &ed_emlrtRSI;
    for (loop_ub = 0; loop_ub < m; loop_ub++) {
      c_st.site = &fd_emlrtRSI;
      if ((1 <= loop_ub + 1) && (loop_ub + 1 > 2147483646)) {
        d_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      for (b_i = 0; b_i <= loop_ub; b_i++) {
        R->data[b_i + R->size[0] * loop_ub] = b_A->data[b_i + b_A->size[0] *
          loop_ub];
      }

      c_st.site = &gd_emlrtRSI;
      if ((loop_ub + 2 <= m) && (m > 2147483646)) {
        d_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      if (loop_ub + 2 <= m) {
        R->data[R->size[0] * loop_ub + 1] = 0.0F;
      }
    }

    i = b_A->size[0] + 1;
    for (loop_ub = i; loop_ub < 3; loop_ub++) {
      c_st.site = &hd_emlrtRSI;
      if ((1 <= m) && (m > 2147483646)) {
        d_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      for (b_i = 0; b_i < m; b_i++) {
        R->data[b_i + R->size[0] * (loop_ub - 1)] = b_A->data[b_i + b_A->size[0]
          * (loop_ub - 1)];
      }
    }

    c_st.site = &id_emlrtRSI;
    i = c_A->size[0] * c_A->size[1];
    c_A->size[0] = b_A->size[0];
    c_A->size[1] = 2;
    emxEnsureCapacity_real32_T(&c_st, c_A, i, &td_emlrtRTEI);
    loop_ub = b_A->size[0] * b_A->size[1];
    for (i = 0; i < loop_ub; i++) {
      c_A->data[i] = b_A->data[i];
    }

    d_st.site = &xd_emlrtRSI;
    e_st.site = &yd_emlrtRSI;
    if (c_A->size[0] != 0) {
      info_t = LAPACKE_sorgqr(102, (ptrdiff_t)b_A->size[0], (ptrdiff_t)b_A->
        size[0], (ptrdiff_t)b_A->size[0], &c_A->data[0], (ptrdiff_t)b_A->size[0],
        &tau_data[0]);
      loop_ub = (int32_T)info_t;
      f_st.site = &de_emlrtRSI;
      if (loop_ub != 0) {
        p = true;
        b_p = false;
        if (loop_ub == -7) {
          b_p = true;
        } else {
          if (loop_ub == -5) {
            b_p = true;
          }
        }

        if (!b_p) {
          if (loop_ub == -1010) {
            emlrtErrorWithMessageIdR2018a(&f_st, &m_emlrtRTEI, "MATLAB:nomem",
              "MATLAB:nomem", 0);
          } else {
            emlrtErrorWithMessageIdR2018a(&f_st, &l_emlrtRTEI,
              "Coder:toolbox:LAPACKCallErrorInfo",
              "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, fname, 12, loop_ub);
          }
        }
      } else {
        p = false;
      }

      if (p) {
        i = c_A->size[0] * c_A->size[1];
        c_A->size[1] = 2;
        emxEnsureCapacity_real32_T(&e_st, c_A, i, &td_emlrtRTEI);
        loop_ub = c_A->size[0];
        for (i = 0; i < loop_ub; i++) {
          c_A->data[i] = rtNaNF;
        }

        loop_ub = c_A->size[0];
        for (i = 0; i < loop_ub; i++) {
          c_A->data[i + c_A->size[0]] = rtNaNF;
        }
      }
    }

    i = Q->size[0] * Q->size[1];
    Q->size[0] = static_cast<int8_T>(b_A->size[0]);
    Q->size[1] = static_cast<int8_T>(b_A->size[0]);
    emxEnsureCapacity_real32_T(&b_st, Q, i, &sd_emlrtRTEI);
    c_st.site = &jd_emlrtRSI;
    for (loop_ub = 0; loop_ub < m; loop_ub++) {
      c_st.site = &kd_emlrtRSI;
      if ((1 <= m) && (m > 2147483646)) {
        d_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      for (b_i = 0; b_i < m; b_i++) {
        Q->data[b_i + Q->size[0] * loop_ub] = c_A->data[b_i + c_A->size[0] *
          loop_ub];
      }
    }
  }

  emxFree_real32_T(&b_A);
  emxFree_real32_T(&c_A);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (qr.cpp) */
