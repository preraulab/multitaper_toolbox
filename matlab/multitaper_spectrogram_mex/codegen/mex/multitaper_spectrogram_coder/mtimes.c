/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.c
 *
 * Code generation for function 'mtimes'
 *
 */

/* Include files */
#include "mtimes.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"
#include "blas.h"
#include <stddef.h>

/* Function Definitions */
real32_T mtimes(const emxArray_real32_T *A, const emxArray_real32_T *B)
{
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  ptrdiff_t n_t;
  real32_T C;
  if (A->size[0] < 1) {
    C = 0.0F;
  } else {
    n_t = (ptrdiff_t)A->size[0];
    incx_t = (ptrdiff_t)1;
    incy_t = (ptrdiff_t)1;
    C = sdot(&n_t, &A->data[0], &incx_t, &B->data[0], &incy_t);
  }

  return C;
}

/* End of code generation (mtimes.c) */
