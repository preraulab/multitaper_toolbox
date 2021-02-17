/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eml_mtimes_helper.h
 *
 * Code generation for function 'eml_mtimes_helper'
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
void b_dynamic_size_checks(const emlrtStack *sp, const emxArray_real32_T *a,
  const emxArray_real_T *b, int32_T innerDimA, int32_T innerDimB);
void dynamic_size_checks(const emlrtStack *sp, const emxArray_real32_T *a, const
  emxArray_real32_T *b, int32_T innerDimA, int32_T innerDimB);

/* End of code generation (eml_mtimes_helper.h) */
