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
static emlrtRSInfo
    cd_emlrtRSI =
        {
            34,       /* lineNo */
            "repmat", /* fcnName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/"
            "repmat.m" /* pathName */
};

static emlrtRSInfo
    dd_emlrtRSI =
        {
            70,       /* lineNo */
            "repmat", /* fcnName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/"
            "repmat.m" /* pathName */
};

static emlrtRSInfo
    ed_emlrtRSI =
        {
            77,       /* lineNo */
            "repmat", /* fcnName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/"
            "repmat.m" /* pathName */
};

static emlrtRTEInfo h_emlrtRTEI = {
    58,                   /* lineNo */
    23,                   /* colNo */
    "assertValidSizeArg", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
    "assertValidSizeArg.m" /* pName */
};

static emlrtRTEInfo
    cd_emlrtRTEI =
        {
            65,       /* lineNo */
            28,       /* colNo */
            "repmat", /* fName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/"
            "repmat.m" /* pName */
};

/* Function Definitions */
void repmat(const emlrtStack *sp, const emxArray_real32_T *a, real_T varargin_2,
            emxArray_real32_T *b)
{
  emlrtStack b_st;
  emlrtStack st;
  int32_T i;
  int32_T i1;
  int32_T ibtile;
  int32_T jtilecol;
  int32_T k;
  const real32_T *a_data;
  real32_T *b_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  a_data = a->data;
  st.site = &cd_emlrtRSI;
  if ((varargin_2 != varargin_2) || muDoubleScalarIsInf(varargin_2)) {
    emlrtErrorWithMessageIdR2018a(
        &st, &h_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  i = a->size[0];
  ibtile = b->size[0] * b->size[1];
  b->size[0] = a->size[0];
  i1 = (int32_T)varargin_2;
  b->size[1] = (int32_T)varargin_2;
  emxEnsureCapacity_real32_T(sp, b, ibtile, &cd_emlrtRTEI);
  b_data = b->data;
  st.site = &dd_emlrtRSI;
  if ((int32_T)varargin_2 > 2147483646) {
    b_st.site = &ab_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }
  for (jtilecol = 0; jtilecol < i1; jtilecol++) {
    ibtile = jtilecol * i;
    st.site = &ed_emlrtRSI;
    if (i > 2147483646) {
      b_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }
    for (k = 0; k < i; k++) {
      b_data[ibtile + k] = a_data[k];
    }
  }
}

/* End of code generation (repmat.c) */
