/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * indexShapeCheck.c
 *
 * Code generation for function 'indexShapeCheck'
 *
 */

/* Include files */
#include "indexShapeCheck.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo db_emlrtRSI = { 33, /* lineNo */
  "indexShapeCheck",                   /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/indexShapeCheck.m"/* pathName */
};

static emlrtRTEInfo f_emlrtRTEI = { 121,/* lineNo */
  5,                                   /* colNo */
  "errOrWarnIf",                       /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/indexShapeCheck.m"/* pName */
};

/* Function Definitions */
void indexShapeCheck(const emlrtStack *sp, const int32_T matrixSize[2], const
                     int32_T indexSize[2])
{
  boolean_T nonSingletonDimFound;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  nonSingletonDimFound = (matrixSize[0] != 1);
  if (matrixSize[1] != 1) {
    if (nonSingletonDimFound) {
      nonSingletonDimFound = false;
    } else {
      nonSingletonDimFound = true;
    }
  }

  if (nonSingletonDimFound) {
    nonSingletonDimFound = (indexSize[0] != 1);
    if (indexSize[1] != 1) {
      if (nonSingletonDimFound) {
        nonSingletonDimFound = false;
      } else {
        nonSingletonDimFound = true;
      }
    }

    if (nonSingletonDimFound) {
      if ((matrixSize[0] != 1) || ((matrixSize[1] == 1) != (indexSize[1] == 1)))
      {
        nonSingletonDimFound = true;
      } else {
        nonSingletonDimFound = false;
      }

      st.site = &db_emlrtRSI;
      if (nonSingletonDimFound) {
        emlrtErrorWithMessageIdR2018a(&st, &f_emlrtRTEI,
          "Coder:FE:PotentialMatrixMatrix_MV",
          "Coder:FE:PotentialMatrixMatrix_MV", 0);
      }
    }
  }
}

/* End of code generation (indexShapeCheck.c) */
