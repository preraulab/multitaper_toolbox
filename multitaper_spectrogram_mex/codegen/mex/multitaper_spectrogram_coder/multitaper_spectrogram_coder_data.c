/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder_data.c
 *
 * Code generation for function 'multitaper_spectrogram_coder_data'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
omp_lock_t emlrtLockGlobal;
omp_nest_lock_t emlrtNestLockGlobal;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131595U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "multitaper_spectrogram_coder",      /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo ab_emlrtRSI = { 21,        /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

emlrtRSInfo gb_emlrtRSI = { 45,        /* lineNo */
  "mpower",                            /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/mpower.m"/* pathName */
};

emlrtRSInfo hb_emlrtRSI = { 70,        /* lineNo */
  "power",                             /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/power.m"/* pathName */
};

emlrtRSInfo jb_emlrtRSI = { 99,        /* lineNo */
  "sumprod",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pathName */
};

emlrtRSInfo ob_emlrtRSI = { 143,       /* lineNo */
  "allOrAny",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/allOrAny.m"/* pathName */
};

emlrtRSInfo xc_emlrtRSI = { 66,        /* lineNo */
  "applyBinaryScalarFunction",         /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

emlrtRSInfo yc_emlrtRSI = { 188,       /* lineNo */
  "flatIter",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

emlrtRSInfo td_emlrtRSI = { 48,        /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pathName */
};

emlrtRSInfo te_emlrtRSI = { 69,        /* lineNo */
  "combineVectorElements",             /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

emlrtRSInfo ue_emlrtRSI = { 87,        /* lineNo */
  "blockedSummation",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRSInfo ve_emlrtRSI = { 149,       /* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRSInfo we_emlrtRSI = { 168,       /* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRTEInfo p_emlrtRTEI = { 123,      /* lineNo */
  23,                                  /* colNo */
  "dynamic_size_checks",               /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pName */
};

emlrtRTEInfo q_emlrtRTEI = { 118,      /* lineNo */
  23,                                  /* colNo */
  "dynamic_size_checks",               /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pName */
};

emlrtRTEInfo id_emlrtRTEI = { 46,      /* lineNo */
  6,                                   /* colNo */
  "applyBinaryScalarFunction",         /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pName */
};

emlrtRTEInfo jd_emlrtRTEI = { 122,     /* lineNo */
  24,                                  /* colNo */
  "blockedSummation",                  /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
};

/* End of code generation (multitaper_spectrogram_coder_data.c) */
