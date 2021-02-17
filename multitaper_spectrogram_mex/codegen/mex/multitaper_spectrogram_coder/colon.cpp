/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * colon.cpp
 *
 * Code generation for function 'colon'
 *
 */

/* Include files */
#include "colon.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo x_emlrtRSI = { 288, /* lineNo */
  "eml_float_colon",                   /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRTEInfo e_emlrtRTEI = { 405,/* lineNo */
  15,                                  /* colNo */
  "assert_pmaxsize",                   /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pName */
};

static emlrtRTEInfo yc_emlrtRTEI = { 279,/* lineNo */
  14,                                  /* colNo */
  "colon",                             /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pName */
};

/* Function Definitions */
void b_eml_float_colon(const emlrtStack *sp, real_T d, real_T b, emxArray_real_T
  *y)
{
  real_T ndbl;
  real_T apnd;
  real_T cdiff;
  int32_T n;
  int32_T nm1d2;
  int32_T k;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  ndbl = muDoubleScalarFloor(b / d + 0.5);
  apnd = ndbl * d;
  if (d > 0.0) {
    cdiff = apnd - b;
  } else {
    cdiff = b - apnd;
  }

  if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarAbs(b))
  {
    ndbl++;
    apnd = b;
  } else if (cdiff > 0.0) {
    apnd = (ndbl - 1.0) * d;
  } else {
    ndbl++;
  }

  if (ndbl >= 0.0) {
    n = static_cast<int32_T>(ndbl);
  } else {
    n = 0;
  }

  st.site = &x_emlrtRSI;
  if (ndbl > 2.147483647E+9) {
    emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI, "Coder:MATLAB:pmaxsize",
      "Coder:MATLAB:pmaxsize", 0);
  }

  nm1d2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity_real_T(sp, y, nm1d2, &yc_emlrtRTEI);
  if (n > 0) {
    y->data[0] = 0.0;
    if (n > 1) {
      y->data[n - 1] = apnd;
      nm1d2 = (n - 1) / 2;
      for (k = 0; k <= nm1d2 - 2; k++) {
        ndbl = (static_cast<real_T>(k) + 1.0) * d;
        y->data[k + 1] = ndbl;
        y->data[(n - k) - 2] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        y->data[nm1d2] = apnd / 2.0;
      } else {
        ndbl = static_cast<real_T>(nm1d2) * d;
        y->data[nm1d2] = ndbl;
        y->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }
}

void eml_float_colon(const emlrtStack *sp, real_T d, real_T b, emxArray_real_T
                     *y)
{
  real_T ndbl;
  real_T apnd;
  real_T cdiff;
  int32_T n;
  int32_T nm1d2;
  int32_T k;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  ndbl = muDoubleScalarFloor((b - 1.0) / d + 0.5);
  apnd = ndbl * d + 1.0;
  if (d > 0.0) {
    cdiff = apnd - b;
  } else {
    cdiff = b - apnd;
  }

  if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax(1.0,
       muDoubleScalarAbs(b))) {
    ndbl++;
    apnd = b;
  } else if (cdiff > 0.0) {
    apnd = (ndbl - 1.0) * d + 1.0;
  } else {
    ndbl++;
  }

  if (ndbl >= 0.0) {
    n = static_cast<int32_T>(ndbl);
  } else {
    n = 0;
  }

  st.site = &x_emlrtRSI;
  if (ndbl > 2.147483647E+9) {
    emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI, "Coder:MATLAB:pmaxsize",
      "Coder:MATLAB:pmaxsize", 0);
  }

  nm1d2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity_real_T(sp, y, nm1d2, &yc_emlrtRTEI);
  if (n > 0) {
    y->data[0] = 1.0;
    if (n > 1) {
      y->data[n - 1] = apnd;
      nm1d2 = (n - 1) / 2;
      for (k = 0; k <= nm1d2 - 2; k++) {
        ndbl = (static_cast<real_T>(k) + 1.0) * d;
        y->data[k + 1] = ndbl + 1.0;
        y->data[(n - k) - 2] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        y->data[nm1d2] = (apnd + 1.0) / 2.0;
      } else {
        ndbl = static_cast<real_T>(nm1d2) * d;
        y->data[nm1d2] = ndbl + 1.0;
        y->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }
}

/* End of code generation (colon.cpp) */
