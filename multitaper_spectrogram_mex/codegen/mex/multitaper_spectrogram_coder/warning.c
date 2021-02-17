/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * warning.c
 *
 * Code generation for function 'warning'
 *
 */

/* Include files */
#include "warning.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo emlrtRTEI = { 14,  /* lineNo */
  15,                                  /* colNo */
  "warning",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/shared/coder/coder/lib/+coder/+internal/warning.m"/* pName */
};

static emlrtMCInfo emlrtMCI = { 14,    /* lineNo */
  25,                                  /* colNo */
  "warning",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/shared/coder/coder/lib/+coder/+internal/warning.m"/* pName */
};

static emlrtMCInfo b_emlrtMCI = { 14,  /* lineNo */
  9,                                   /* colNo */
  "warning",                           /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/shared/coder/coder/lib/+coder/+internal/warning.m"/* pName */
};

static emlrtRSInfo yf_emlrtRSI = { 14, /* lineNo */
  "warning",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/shared/coder/coder/lib/+coder/+internal/warning.m"/* pathName */
};

/* Function Declarations */
static void b_feval(const emlrtStack *sp, const mxArray *b, const mxArray *c,
                    emlrtMCInfo *location);
static const mxArray *feval(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location);

/* Function Definitions */
static void b_feval(const emlrtStack *sp, const mxArray *b, const mxArray *c,
                    emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  pArrays[0] = b;
  pArrays[1] = c;
  emlrtCallMATLABR2012b(sp, 0, NULL, 2, pArrays, "feval", true, location);
}

static const mxArray *feval(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(sp, 1, &m, 2, pArrays, "feval", true, location);
}

void warning(const emlrtStack *sp)
{
  static const int32_T iv[2] = { 1, 7 };

  static const int32_T iv1[2] = { 1, 7 };

  static const int32_T iv2[2] = { 1, 28 };

  static const char_T msgID[28] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'd', 'e',
    't', 'r', 'e', 'n', 'd', ':', 'P', 'o', 'l', 'y', 'N', 'o', 't', 'U', 'n',
    'i', 'q', 'u', 'e' };

  static const char_T b_u[7] = { 'm', 'e', 's', 's', 'a', 'g', 'e' };

  static const char_T u[7] = { 'w', 'a', 'r', 'n', 'i', 'n', 'g' };

  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  const mxArray *y;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtAssertMATLABThread(sp, &emlrtRTEI);
  y = NULL;
  m = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a(sp, 7, m, &u[0]);
  emlrtAssign(&y, m);
  b_y = NULL;
  m = emlrtCreateCharArray(2, &iv1[0]);
  emlrtInitCharArrayR2013a(sp, 7, m, &b_u[0]);
  emlrtAssign(&b_y, m);
  c_y = NULL;
  m = emlrtCreateCharArray(2, &iv2[0]);
  emlrtInitCharArrayR2013a(sp, 28, m, &msgID[0]);
  emlrtAssign(&c_y, m);
  st.site = &yf_emlrtRSI;
  b_feval(&st, y, feval(&st, b_y, c_y, &emlrtMCI), &b_emlrtMCI);
}

/* End of code generation (warning.c) */
