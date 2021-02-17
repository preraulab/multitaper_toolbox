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
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo ge_emlrtRSI = { 28, /* lineNo */
  "repmat",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/repmat.m"/* pathName */
};

static emlrtRSInfo he_emlrtRSI = { 64, /* lineNo */
  "repmat",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/repmat.m"/* pathName */
};

static emlrtRSInfo ie_emlrtRSI = { 71, /* lineNo */
  "repmat",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/repmat.m"/* pathName */
};

static emlrtRTEInfo l_emlrtRTEI = { 58,/* lineNo */
  23,                                  /* colNo */
  "assertValidSizeArg",                /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/assertValidSizeArg.m"/* pName */
};

static emlrtRTEInfo fd_emlrtRTEI = { 59,/* lineNo */
  28,                                  /* colNo */
  "repmat",                            /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/repmat.m"/* pName */
};

/* Function Definitions */
void repmat(const emlrtStack *sp, const emxArray_real32_T *a, real_T varargin_2,
            emxArray_real32_T *b)
{
  emlrtStack b_st;
  emlrtStack st;
  int32_T ibtile;
  int32_T jtilecol;
  int32_T k;
  int32_T nrows;
  int32_T ntilecols;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &ge_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  if ((varargin_2 != varargin_2) || muDoubleScalarIsInf(varargin_2)) {
    emlrtErrorWithMessageIdR2018a(&st, &l_emlrtRTEI,
      "Coder:MATLAB:NonIntegerInput", "Coder:MATLAB:NonIntegerInput", 4, 12,
      MIN_int32_T, 12, MAX_int32_T);
  }

  nrows = b->size[0] * b->size[1];
  b->size[0] = a->size[0];
  b->size[1] = (int32_T)varargin_2;
  emxEnsureCapacity_real32_T(sp, b, nrows, &fd_emlrtRTEI);
  nrows = a->size[0];
  ntilecols = (int32_T)varargin_2;
  st.site = &he_emlrtRSI;
  if ((1 <= (int32_T)varargin_2) && ((int32_T)varargin_2 > 2147483646)) {
    b_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }

  for (jtilecol = 0; jtilecol < ntilecols; jtilecol++) {
    ibtile = jtilecol * nrows;
    st.site = &ie_emlrtRSI;
    if ((1 <= nrows) && (nrows > 2147483646)) {
      b_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (k = 0; k < nrows; k++) {
      b->data[ibtile + k] = a->data[k];
    }
  }
}

/* End of code generation (repmat.c) */
