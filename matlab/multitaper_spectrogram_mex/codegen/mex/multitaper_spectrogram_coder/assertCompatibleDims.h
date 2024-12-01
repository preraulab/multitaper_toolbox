/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * assertCompatibleDims.h
 *
 * Code generation for function 'assertCompatibleDims'
 *
 */

#pragma once

/* Include files */
#include "multitaper_spectrogram_coder_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void assertCompatibleDims(const emlrtStack *sp, const emxArray_real32_T *x,
                          const emxArray_real32_T *y);

void b_assertCompatibleDims(const emlrtStack *sp, const emxArray_real32_T *x,
                            const emxArray_real32_T *y);

/* End of code generation (assertCompatibleDims.h) */
