/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * blockedSummation.h
 *
 * Code generation for function 'blockedSummation'
 *
 */

#pragma once

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "omp.h"
#include "multitaper_spectrogram_coder_mex_types.h"

/* Function Declarations */
void blockedSummation(const emlrtStack *sp, const emxArray_real32_T *x, int32_T
                      vlen, emxArray_real32_T *y);

/* End of code generation (blockedSummation.h) */
