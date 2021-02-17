/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.cpp
 *
 * Code generation for function 'mtimes'
 *
 */

/* Include files */
#include "mtimes.h"
#include "blas.h"
#include "multitaper_spectrogram_coder.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo ce_emlrtRTEI = { 166,/* lineNo */
  9,                                   /* colNo */
  "mtimes",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pName */
};

/* Function Definitions */
real32_T mtimes(const emlrtStack *sp, const emxArray_real32_T *A, const
                emxArray_real32_T *B)
{
  real32_T C;
  emxArray_real32_T *At;
  int32_T i;
  int32_T loop_ub;
  ptrdiff_t n_t;
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real32_T(sp, &At, 2, &ce_emlrtRTEI, true);
  i = At->size[0] * At->size[1];
  At->size[0] = 1;
  At->size[1] = A->size[0];
  emxEnsureCapacity_real32_T(sp, At, i, &ce_emlrtRTEI);
  loop_ub = A->size[0];
  for (i = 0; i < loop_ub; i++) {
    At->data[i] = A->data[i];
  }

  if (A->size[0] < 1) {
    C = 0.0F;
  } else {
    n_t = (ptrdiff_t)A->size[0];
    incx_t = (ptrdiff_t)1;
    incy_t = (ptrdiff_t)1;
    C = sdot(&n_t, &At->data[0], &incx_t, &B->data[0], &incy_t);
  }

  emxFree_real32_T(&At);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
  return C;
}

/* End of code generation (mtimes.cpp) */
