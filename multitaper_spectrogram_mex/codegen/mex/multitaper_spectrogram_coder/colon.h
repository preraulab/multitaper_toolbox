/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * colon.h
 *
 * Code generation for function 'colon'
 *
 */

#pragma once

/* Include files */
#include "multitaper_spectrogram_coder_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include "omp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void b_eml_float_colon(const emlrtStack *sp, real_T d, real_T b, emxArray_real_T
  *y);
void eml_float_colon(const emlrtStack *sp, real_T d, real_T b, emxArray_real_T
                     *y);
void eml_integer_colon_dispatcher(const emlrtStack *sp, int32_T b,
  emxArray_int32_T *y);

/* End of code generation (colon.h) */
