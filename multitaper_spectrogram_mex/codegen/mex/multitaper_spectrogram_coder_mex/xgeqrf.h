/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgeqrf.h
 *
 * Code generation for function 'xgeqrf'
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
void xgeqrf(const emlrtStack *sp, emxArray_real32_T *A, real32_T tau_data[],
            int32_T tau_size[1]);

/* End of code generation (xgeqrf.h) */
