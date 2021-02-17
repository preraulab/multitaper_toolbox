/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * trisolve.c
 *
 * Code generation for function 'trisolve'
 *
 */

/* Include files */
#include "trisolve.h"
#include "rt_nonfinite.h"
#include "blas.h"
#include <stddef.h>

/* Function Definitions */
void trisolve(const real32_T A_data[], const int32_T A_size[2], real32_T B[2])
{
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  real32_T alpha1;
  char_T DIAGA1;
  char_T SIDE1;
  char_T TRANSA1;
  char_T UPLO1;
  if (A_size[0] >= 1) {
    alpha1 = 1.0F;
    DIAGA1 = 'N';
    TRANSA1 = 'N';
    UPLO1 = 'U';
    SIDE1 = 'L';
    m_t = (ptrdiff_t)A_size[0];
    n_t = (ptrdiff_t)1;
    lda_t = (ptrdiff_t)A_size[0];
    ldb_t = (ptrdiff_t)2;
    strsm(&SIDE1, &UPLO1, &TRANSA1, &DIAGA1, &m_t, &n_t, &alpha1, &A_data[0],
          &lda_t, &B[0], &ldb_t);
  }
}

/* End of code generation (trisolve.c) */
