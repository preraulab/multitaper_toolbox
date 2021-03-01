/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * trisolve.cpp
 *
 * Code generation for function 'trisolve'
 *
 */

/* Include files */
#include "trisolve.h"
#include "blas.h"
#include "multitaper_spectrogram_coder.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void trisolve(const emxArray_real32_T *A, real32_T B[2])
{
  int32_T mA;
  real32_T alpha1;
  char_T DIAGA1;
  char_T TRANSA1;
  char_T UPLO1;
  char_T SIDE1;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  mA = A->size[0];
  mA = muIntScalarMin_sint32(mA, 2);
  if (mA >= 1) {
    alpha1 = 1.0F;
    DIAGA1 = 'N';
    TRANSA1 = 'N';
    UPLO1 = 'U';
    SIDE1 = 'L';
    m_t = (ptrdiff_t)mA;
    n_t = (ptrdiff_t)1;
    lda_t = (ptrdiff_t)A->size[0];
    ldb_t = (ptrdiff_t)2;
    strsm(&SIDE1, &UPLO1, &TRANSA1, &DIAGA1, &m_t, &n_t, &alpha1, &A->data[0],
          &lda_t, &B[0], &ldb_t);
  }
}

/* End of code generation (trisolve.cpp) */
