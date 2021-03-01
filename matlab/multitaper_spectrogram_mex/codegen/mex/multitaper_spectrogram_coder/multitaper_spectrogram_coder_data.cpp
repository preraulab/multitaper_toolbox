/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder_data.cpp
 *
 * Code generation for function 'multitaper_spectrogram_coder_data'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder.h"
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
  "multitaper_spectrogram_coder",      /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo y_emlrtRSI = { 306,        /* lineNo */
  "eml_float_colon",                   /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

emlrtRSInfo ab_emlrtRSI = { 21,        /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

emlrtRSInfo bb_emlrtRSI = { 11,        /* lineNo */
  "nextpow2",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elfun/nextpow2.m"/* pathName */
};

emlrtRSInfo cb_emlrtRSI = { 17,        /* lineNo */
  "applyScalarFunctionInPlace",        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/applyScalarFunctionInPlace.m"/* pathName */
};

emlrtRSInfo db_emlrtRSI = { 11,        /* lineNo */
  "nextpow2",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+scalar/nextpow2.m"/* pathName */
};

emlrtRSInfo eb_emlrtRSI = { 23,        /* lineNo */
  "fnextpow2",                         /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+scalar/nextpow2.m"/* pathName */
};

emlrtRSInfo fb_emlrtRSI = { 17,        /* lineNo */
  "log2",                              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+scalar/log2.m"/* pathName */
};

emlrtRSInfo gb_emlrtRSI = { 45,        /* lineNo */
  "mpower",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/mpower.m"/* pathName */
};

emlrtRSInfo hb_emlrtRSI = { 55,        /* lineNo */
  "power",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/power.m"/* pathName */
};

emlrtRSInfo jb_emlrtRSI = { 96,        /* lineNo */
  "sumprod",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pathName */
};

emlrtRSInfo nb_emlrtRSI = { 143,       /* lineNo */
  "allOrAny",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/allOrAny.m"/* pathName */
};

emlrtRSInfo wc_emlrtRSI = { 66,        /* lineNo */
  "applyBinaryScalarFunction",         /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

emlrtRSInfo xc_emlrtRSI = { 188,       /* lineNo */
  "flatIter",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

emlrtRSInfo ud_emlrtRSI = { 9,         /* lineNo */
  "int",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/int.m"/* pathName */
};

emlrtRSInfo vd_emlrtRSI = { 8,         /* lineNo */
  "majority",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/majority.m"/* pathName */
};

emlrtRSInfo wd_emlrtRSI = { 31,        /* lineNo */
  "infocheck",                         /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pathName */
};

emlrtRSInfo ae_emlrtRSI = { 38,        /* lineNo */
  "ceval_xorgqr",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

emlrtRSInfo be_emlrtRSI = { 46,        /* lineNo */
  "ceval_xorgqr",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

emlrtRSInfo ce_emlrtRSI = { 51,        /* lineNo */
  "ceval_xorgqr",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xorgqr.m"/* pathName */
};

emlrtRSInfo ee_emlrtRSI = { 102,       /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pathName */
};

emlrtRSInfo fe_emlrtRSI = { 56,        /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pathName */
};

emlrtRSInfo me_emlrtRSI = { 99,        /* lineNo */
  "trisolve",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/trisolve.m"/* pathName */
};

emlrtRSInfo ne_emlrtRSI = { 122,       /* lineNo */
  "trisolveBLAS",                      /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/trisolve.m"/* pathName */
};

emlrtRSInfo oe_emlrtRSI = { 77,        /* lineNo */
  "xtrsm",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/xtrsm.m"/* pathName */
};

emlrtRSInfo pe_emlrtRSI = { 76,        /* lineNo */
  "xtrsm",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/xtrsm.m"/* pathName */
};

emlrtRSInfo we_emlrtRSI = { 18,        /* lineNo */
  "fftw",                              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/fftw.m"/* pathName */
};

emlrtRSInfo xe_emlrtRSI = { 28,        /* lineNo */
  "FFTWApi/fft1d",                     /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+fftw/FFTWApi.m"/* pathName */
};

emlrtRSInfo af_emlrtRSI = { 15,        /* lineNo */
  "MATLABFFTWCallback/fft1d",          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+fftw/MATLABFFTWCallback.m"/* pathName */
};

emlrtRSInfo bf_emlrtRSI = { 17,        /* lineNo */
  "MATLABFFTWCallback/fft1d",          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+fftw/MATLABFFTWCallback.m"/* pathName */
};

emlrtRSInfo df_emlrtRSI = { 174,       /* lineNo */
  "mtimes",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pathName */
};

emlrtRSInfo ef_emlrtRSI = { 32,        /* lineNo */
  "xdotu",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/xdotu.m"/* pathName */
};

emlrtRSInfo ff_emlrtRSI = { 49,        /* lineNo */
  "xdot",                              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/xdot.m"/* pathName */
};

emlrtRSInfo hf_emlrtRSI = { 62,        /* lineNo */
  "combineVectorElements",             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

emlrtRSInfo if_emlrtRSI = { 81,        /* lineNo */
  "blockedSummation",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRSInfo jf_emlrtRSI = { 142,       /* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRSInfo kf_emlrtRSI = { 161,       /* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRSInfo lf_emlrtRSI = { 20,        /* lineNo */
  "sum",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

emlrtRSInfo nf_emlrtRSI = { 159,       /* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRSInfo of_emlrtRSI = { 173,       /* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRSInfo pf_emlrtRSI = { 194,       /* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pathName */
};

emlrtRSInfo yf_emlrtRSI = { 224,       /* lineNo */
  "do_vectors",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

emlrtRSInfo ag_emlrtRSI = { 227,       /* lineNo */
  "do_vectors",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

emlrtRSInfo bg_emlrtRSI = { 258,       /* lineNo */
  "do_vectors",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

emlrtRSInfo cg_emlrtRSI = { 259,       /* lineNo */
  "do_vectors",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

emlrtRSInfo dg_emlrtRSI = { 261,       /* lineNo */
  "do_vectors",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

emlrtRSInfo eg_emlrtRSI = { 348,       /* lineNo */
  "do_vectors",                        /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

emlrtRSInfo fg_emlrtRSI = { 71,        /* lineNo */
  "issorted",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/issorted.m"/* pathName */
};

emlrtRSInfo gg_emlrtRSI = { 93,        /* lineNo */
  "looper",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/issorted.m"/* pathName */
};

emlrtRSInfo hg_emlrtRSI = { 459,       /* lineNo */
  "skip_to_last_equal_value",          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

emlrtRSInfo ig_emlrtRSI = { 40,        /* lineNo */
  "safeEq",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/safeEq.m"/* pathName */
};

emlrtRSInfo jg_emlrtRSI = { 46,        /* lineNo */
  "eps",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elmat/eps.m"/* pathName */
};

emlrtRTEInfo g_emlrtRTEI = { 47,       /* lineNo */
  19,                                  /* colNo */
  "allOrAny",                          /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/allOrAny.m"/* pName */
};

emlrtRTEInfo h_emlrtRTEI = { 34,       /* lineNo */
  23,                                  /* colNo */
  "detrend",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/detrend.m"/* pName */
};

emlrtRTEInfo j_emlrtRTEI = { 153,      /* lineNo */
  23,                                  /* colNo */
  "dynamic_size_checks",               /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pName */
};

emlrtRTEInfo k_emlrtRTEI = { 158,      /* lineNo */
  23,                                  /* colNo */
  "dynamic_size_checks",               /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pName */
};

emlrtRTEInfo l_emlrtRTEI = { 48,       /* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

emlrtRTEInfo m_emlrtRTEI = { 45,       /* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

emlrtRTEInfo de_emlrtRTEI = { 81,      /* lineNo */
  13,                                  /* colNo */
  "blockedSummation",                  /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/private/blockedSummation.m"/* pName */
};

/* End of code generation (multitaper_spectrogram_coder_data.cpp) */
