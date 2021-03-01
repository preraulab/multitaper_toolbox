/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgeqrf.cpp
 *
 * Code generation for function 'xgeqrf'
 *
 */

/* Include files */
#include "xgeqrf.h"
#include "eml_int_forloop_overflow_check.h"
#include "lapacke.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_data.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo ld_emlrtRSI = { 27, /* lineNo */
  "xgeqrf",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo md_emlrtRSI = { 52, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo nd_emlrtRSI = { 53, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo od_emlrtRSI = { 79, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo pd_emlrtRSI = { 84, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo qd_emlrtRSI = { 91, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo rd_emlrtRSI = { 94, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo sd_emlrtRSI = { 99, /* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

static emlrtRSInfo td_emlrtRSI = { 102,/* lineNo */
  "ceval_xgeqrf",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgeqrf.m"/* pathName */
};

/* Function Definitions */
void xgeqrf(const emlrtStack *sp, emxArray_real32_T *A, real32_T tau_data[],
            int32_T tau_size[1])
{
  int32_T m;
  int32_T ma;
  int32_T minmana;
  ptrdiff_t info_t;
  int32_T minmn;
  boolean_T p;
  static const char_T fname[14] = { 'L', 'A', 'P', 'A', 'C', 'K', 'E', '_', 's',
    'g', 'e', 'q', 'r', 'f' };

  int32_T i;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  m = A->size[0];
  st.site = &ld_emlrtRSI;
  ma = A->size[0];
  minmana = muIntScalarMin_sint32(ma, 2);
  b_st.site = &md_emlrtRSI;
  b_st.site = &nd_emlrtRSI;
  tau_size[0] = minmana;
  if (A->size[0] == 0) {
    tau_size[0] = minmana;
    if (0 <= minmana - 1) {
      memset(&tau_data[0], 0, minmana * sizeof(real32_T));
    }
  } else {
    b_st.site = &od_emlrtRSI;
    b_st.site = &pd_emlrtRSI;
    info_t = LAPACKE_sgeqrf(102, (ptrdiff_t)A->size[0], (ptrdiff_t)2, &A->data[0],
      (ptrdiff_t)A->size[0], &tau_data[0]);
    minmn = (int32_T)info_t;
    b_st.site = &qd_emlrtRSI;
    if (minmn != 0) {
      p = true;
      if (minmn != -4) {
        if (minmn == -1010) {
          emlrtErrorWithMessageIdR2018a(&b_st, &m_emlrtRTEI, "MATLAB:nomem",
            "MATLAB:nomem", 0);
        } else {
          emlrtErrorWithMessageIdR2018a(&b_st, &l_emlrtRTEI,
            "Coder:toolbox:LAPACKCallErrorInfo",
            "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, fname, 12, minmn);
        }
      }
    } else {
      p = false;
    }

    if (p) {
      b_st.site = &rd_emlrtRSI;
      if ((1 <= m) && (m > 2147483646)) {
        c_st.site = &ab_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      for (i = 0; i < m; i++) {
        A->data[i] = rtNaNF;
      }

      b_st.site = &rd_emlrtRSI;
      for (i = 0; i < m; i++) {
        A->data[ma + i] = rtNaNF;
      }

      minmn = muIntScalarMin_sint32(m, 2) - 1;
      b_st.site = &sd_emlrtRSI;
      for (i = 0; i <= minmn; i++) {
        tau_data[i] = rtNaNF;
      }

      b_st.site = &td_emlrtRSI;
      if (minmn + 2 <= minmana) {
        tau_data[1] = 0.0F;
      }
    }
  }
}

/* End of code generation (xgeqrf.cpp) */
