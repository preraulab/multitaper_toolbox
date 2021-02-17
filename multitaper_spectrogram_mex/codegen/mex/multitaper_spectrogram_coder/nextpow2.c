/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nextpow2.c
 *
 * Code generation for function 'nextpow2'
 *
 */

/* Include files */
#include "nextpow2.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"
#include <math.h>

/* Function Definitions */
real_T nextpow2(real_T n)
{
  real_T f;
  real_T p;
  int32_T eint;
  p = muDoubleScalarAbs(n);
  if ((!muDoubleScalarIsInf(p)) && (!muDoubleScalarIsNaN(p))) {
    f = frexp(p, &eint);
    p = eint;
    if (f == 0.5) {
      p = (real_T)eint - 1.0;
    }
  }

  return p;
}

/* End of code generation (nextpow2.c) */
