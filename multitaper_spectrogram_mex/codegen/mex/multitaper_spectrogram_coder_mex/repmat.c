/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 */

/* Include files */
#include "repmat.h"
#include "eml_int_forloop_overflow_check.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex_emxutil.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo je_emlrtRSI = { 28, /* lineNo */
  "repmat",                            /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/elmat/repmat.m"/* pathName */
};

static emlrtRSInfo ke_emlrtRSI = { 64, /* lineNo */
  "repmat",                            /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/elmat/repmat.m"/* pathName */
};

static emlrtRSInfo le_emlrtRSI = { 71, /* lineNo */
  "repmat",                            /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/elmat/repmat.m"/* pathName */
};

static emlrtDCInfo c_emlrtDCI = { 31,  /* lineNo */
  14,                                  /* colNo */
  "repmat",                            /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/elmat/repmat.m",/* pName */
  4                                    /* checkKind */
};

static emlrtRTEInfo m_emlrtRTEI = { 58,/* lineNo */
  23,                                  /* colNo */
  "assertValidSizeArg",                /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/assertValidSizeArg.m"/* pName */
};

static emlrtRTEInfo n_emlrtRTEI = { 64,/* lineNo */
  15,                                  /* colNo */
  "assertValidSizeArg",                /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/assertValidSizeArg.m"/* pName */
};

static emlrtRTEInfo pc_emlrtRTEI = { 1,/* lineNo */
  14,                                  /* colNo */
  "repmat",                            /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/elmat/repmat.m"/* pName */
};

/* Function Definitions */
void repmat(const emlrtStack *sp, const emxArray_real32_T *a, real_T varargin_2,
            emxArray_real32_T *b)
{
  real_T b_varargin_2;
  int32_T nrows;
  int32_T i;
  int32_T jtilecol;
  int32_T ibtile;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &je_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  if ((varargin_2 != muDoubleScalarFloor(varargin_2)) || muDoubleScalarIsInf
      (varargin_2) || (varargin_2 < -2.147483648E+9) || (varargin_2 >
       2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&st, &m_emlrtRTEI,
      "Coder:MATLAB:NonIntegerInput", "Coder:MATLAB:NonIntegerInput", 4, 12,
      MIN_int32_T, 12, MAX_int32_T);
  }

  if (varargin_2 <= 0.0) {
    b_varargin_2 = 0.0;
  } else {
    b_varargin_2 = varargin_2;
  }

  if (!(b_varargin_2 <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&st, &n_emlrtRTEI, "Coder:MATLAB:pmaxsize",
      "Coder:MATLAB:pmaxsize", 0);
  }

  if (!(varargin_2 >= 0.0)) {
    emlrtNonNegativeCheckR2012b(varargin_2, &c_emlrtDCI, sp);
  }

  nrows = b->size[0] * b->size[1];
  b->size[0] = a->size[0];
  i = (int32_T)varargin_2;
  b->size[1] = i;
  emxEnsureCapacity_real32_T(sp, b, nrows, &pc_emlrtRTEI);
  nrows = a->size[0];
  st.site = &ke_emlrtRSI;
  if ((1 <= i) && (i > 2147483646)) {
    b_st.site = &q_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }

  for (jtilecol = 0; jtilecol < i; jtilecol++) {
    ibtile = jtilecol * nrows;
    st.site = &le_emlrtRSI;
    if ((1 <= nrows) && (nrows > 2147483646)) {
      b_st.site = &q_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (k = 0; k < nrows; k++) {
      b->data[ibtile + k] = a->data[k];
    }
  }
}

/* End of code generation (repmat.c) */
