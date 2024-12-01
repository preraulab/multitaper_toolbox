/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * assertCompatibleDims.c
 *
 * Code generation for function 'assertCompatibleDims'
 *
 */

/* Include files */
#include "assertCompatibleDims.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo n_emlrtRTEI = {
    13,                     /* lineNo */
    27,                     /* colNo */
    "assertCompatibleDims", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/shared/coder/coder/lib/+coder/"
    "+internal/assertCompatibleDims.m" /* pName */
};

/* Function Definitions */
void assertCompatibleDims(const emlrtStack *sp, const emxArray_real32_T *x,
                          const emxArray_real32_T *y)
{
  if (((x->size[0] != 1) && (y->size[0] != 1) && (x->size[0] != y->size[0])) ||
      ((x->size[1] != 1) && (y->size[1] != 1) && (x->size[1] != y->size[1]))) {
    emlrtErrorWithMessageIdR2018a(sp, &n_emlrtRTEI,
                                  "MATLAB:sizeDimensionsMustMatch",
                                  "MATLAB:sizeDimensionsMustMatch", 0);
  }
}

void b_assertCompatibleDims(const emlrtStack *sp, const emxArray_real32_T *x,
                            const emxArray_real32_T *y)
{
  if ((x->size[0] != 1) && (y->size[0] != 1) && (x->size[0] != y->size[0])) {
    emlrtErrorWithMessageIdR2018a(sp, &n_emlrtRTEI,
                                  "MATLAB:sizeDimensionsMustMatch",
                                  "MATLAB:sizeDimensionsMustMatch", 0);
  }
}

/* End of code generation (assertCompatibleDims.c) */
