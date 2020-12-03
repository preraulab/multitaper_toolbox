/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder_mex_data.c
 *
 * Code generation for function 'multitaper_spectrogram_coder_mex_data'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
omp_lock_t emlrtLockGlobal;
omp_nest_lock_t emlrtNestLockGlobal;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131483U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "multitaper_spectrogram_coder_mex",  /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo p_emlrtRSI = { 306,        /* lineNo */
  "eml_float_colon",                   /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

emlrtRSInfo q_emlrtRSI = { 21,         /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

emlrtRSInfo r_emlrtRSI = { 11,         /* lineNo */
  "nextpow2",                          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/elfun/nextpow2.m"/* pathName */
};

emlrtRSInfo s_emlrtRSI = { 17,         /* lineNo */
  "applyScalarFunctionInPlace",        /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/applyScalarFunctionInPlace.m"/* pathName */
};

emlrtRSInfo t_emlrtRSI = { 11,         /* lineNo */
  "nextpow2",                          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+scalar/nextpow2.m"/* pathName */
};

emlrtRSInfo u_emlrtRSI = { 23,         /* lineNo */
  "fnextpow2",                         /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+scalar/nextpow2.m"/* pathName */
};

emlrtRSInfo v_emlrtRSI = { 17,         /* lineNo */
  "log2",                              /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+scalar/log2.m"/* pathName */
};

emlrtRSInfo w_emlrtRSI = { 45,         /* lineNo */
  "mpower",                            /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/ops/mpower.m"/* pathName */
};

emlrtRSInfo x_emlrtRSI = { 55,         /* lineNo */
  "power",                             /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/ops/power.m"/* pathName */
};

emlrtRSInfo oc_emlrtRSI = { 66,        /* lineNo */
  "applyBinaryScalarFunction",         /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

emlrtRSInfo pc_emlrtRSI = { 188,       /* lineNo */
  "flatIter",                          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

emlrtRSInfo md_emlrtRSI = { 9,         /* lineNo */
  "int",                               /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+lapack/int.m"/* pathName */
};

emlrtRSInfo nd_emlrtRSI = { 8,         /* lineNo */
  "majority",                          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+lapack/majority.m"/* pathName */
};

emlrtRSInfo od_emlrtRSI = { 31,        /* lineNo */
  "infocheck",                         /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pathName */
};

emlrtRSInfo rd_emlrtRSI = { 38,        /* lineNo */
  "ceval_xorgqr",                      /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

emlrtRSInfo sd_emlrtRSI = { 46,        /* lineNo */
  "ceval_xorgqr",                      /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

emlrtRSInfo td_emlrtRSI = { 51,        /* lineNo */
  "ceval_xorgqr",                      /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

emlrtRSInfo ee_emlrtRSI = { 99,        /* lineNo */
  "trisolve",                          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/trisolve.m"/* pathName */
};

emlrtRSInfo fe_emlrtRSI = { 122,       /* lineNo */
  "trisolveBLAS",                      /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/trisolve.m"/* pathName */
};

emlrtRSInfo ge_emlrtRSI = { 77,        /* lineNo */
  "xtrsm",                             /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+blas/xtrsm.m"/* pathName */
};

emlrtRSInfo he_emlrtRSI = { 76,        /* lineNo */
  "xtrsm",                             /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+blas/xtrsm.m"/* pathName */
};

emlrtRSInfo pe_emlrtRSI = { 18,        /* lineNo */
  "fftw",                              /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/fftw.m"/* pathName */
};

emlrtRSInfo se_emlrtRSI = { 15,        /* lineNo */
  "MATLABFFTWCallback/fft1d",          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+fftw/MATLABFFTWCallback.m"/* pathName */
};

emlrtRSInfo te_emlrtRSI = { 17,        /* lineNo */
  "MATLABFFTWCallback/fft1d",          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+fftw/MATLABFFTWCallback.m"/* pathName */
};

emlrtRTEInfo b_emlrtRTEI = { 43,       /* lineNo */
  23,                                  /* colNo */
  "sumprod",                           /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pName */
};

emlrtRTEInfo c_emlrtRTEI = { 47,       /* lineNo */
  19,                                  /* colNo */
  "allOrAny",                          /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/allOrAny.m"/* pName */
};

emlrtRTEInfo g_emlrtRTEI = { 34,       /* lineNo */
  23,                                  /* colNo */
  "detrend",                           /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

emlrtRTEInfo k_emlrtRTEI = { 48,       /* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

emlrtRTEInfo l_emlrtRTEI = { 45,       /* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

emlrtDCInfo d_emlrtDCI = { 19,         /* lineNo */
  34,                                  /* colNo */
  "allocFftOutput",                    /* fName */
  "/apps/source/matlab/2019b/toolbox/coder/coder/+coder/+fftw/allocFftOutput.m",/* pName */
  4                                    /* checkKind */
};

/* End of code generation (multitaper_spectrogram_coder_mex_data.c) */
