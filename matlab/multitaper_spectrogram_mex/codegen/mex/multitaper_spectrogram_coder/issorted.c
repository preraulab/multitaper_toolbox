/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * issorted.c
 *
 * Code generation for function 'issorted'
 *
 */

/* Include files */
#include "issorted.h"
#include "multitaper_spectrogram_coder_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
boolean_T issorted(const emxArray_real_T *x)
{
  real_T v_idx_1;
  int32_T k;
  boolean_T exitg1;
  boolean_T y;
  y = true;
  if ((x->size[1] != 0) && (x->size[1] != 1)) {
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= x->size[1] - 2)) {
      v_idx_1 = x->data[k + 1];
      if ((!(x->data[k] <= v_idx_1)) && (!muDoubleScalarIsNaN(v_idx_1))) {
        y = false;
      }

      if (!y) {
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  return y;
}

/* End of code generation (issorted.c) */
