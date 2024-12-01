/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * div.h
 *
 * Code generation for function 'div'
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
void b_rdivide(const emlrtStack *sp, emxArray_real32_T *in1,
               const emxArray_real32_T *in2);

void rdivide(const emlrtStack *sp, emxArray_real32_T *in1,
             const emxArray_real32_T *in2);

/* End of code generation (div.h) */
