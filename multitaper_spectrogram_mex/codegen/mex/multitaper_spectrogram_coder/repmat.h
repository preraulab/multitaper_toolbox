/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.h
 *
 * Code generation for function 'repmat'
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
void repmat(const emlrtStack *sp, const emxArray_real32_T *a, real_T varargin_2,
            emxArray_real32_T *b);

/* End of code generation (repmat.h) */
