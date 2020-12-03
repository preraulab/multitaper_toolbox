/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * detrend.h
 *
 * Code generation for function 'detrend'
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
void b_detrend(const emlrtStack *sp, const emxArray_real32_T *x,
               emxArray_real32_T *y);
void detrend(const emlrtStack *sp, const emxArray_real32_T *x, emxArray_real32_T
             *y);

/* End of code generation (detrend.h) */
