/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder.cpp
 *
 * Code generation for function 'multitaper_spectrogram_coder'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder.h"
#include "all.h"
#include "any1.h"
#include "blockedSummation.h"
#include "colon.h"
#include "combineVectorElements.h"
#include "detrend.h"
#include "eml_setop.h"
#include "fft.h"
#include "find.h"
#include "mean.h"
#include "mtimes.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "mwmathutil.h"
#include "nextpow2.h"
#include "power.h"
#include "repmat.h"
#include "rt_nonfinite.h"
#include "scalexpCompatible.h"
#include "sum.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 55,    /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 61,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 65,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 77,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 92,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 96,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 103, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 105, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 109, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 112, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 115, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 120, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 121, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 125, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 127, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 128, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 134, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 138, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 146, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 147, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 148, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 149, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 103, /* lineNo */
  "colon",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 15, /* lineNo */
  "sum",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

static emlrtRSInfo tf_emlrtRSI = { 41, /* lineNo */
  "find",                              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elmat/find.m"/* pathName */
};

static emlrtRSInfo wf_emlrtRSI = { 19, /* lineNo */
  "setdiff",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/setdiff.m"/* pathName */
};

static emlrtRSInfo xf_emlrtRSI = { 70, /* lineNo */
  "eml_setop",                         /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

static emlrtRSInfo kg_emlrtRSI = { 27, /* lineNo */
  "cat",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/cat.m"/* pathName */
};

static emlrtRSInfo lg_emlrtRSI = { 102,/* lineNo */
  "cat_impl",                          /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/cat.m"/* pathName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  149,                                 /* lineNo */
  91,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  149,                                 /* lineNo */
  34,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 89,    /* lineNo */
  25,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  89,                                  /* lineNo */
  25,                                  /* colNo */
  "data",                              /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  89,                                  /* lineNo */
  25,                                  /* colNo */
  "window_start",                      /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 149, /* lineNo */
  65,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  149,                                 /* lineNo */
  65,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  69,                                  /* lineNo */
  10,                                  /* colNo */
  "sfreqs",                            /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = { 77,  /* lineNo */
  24,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  4                                    /* checkKind */
};

static emlrtECInfo emlrtECI = { 2,     /* nDims */
  68,                                  /* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtECInfo b_emlrtECI = { 2,   /* nDims */
  109,                                 /* lineNo */
  20,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtECInfo c_emlrtECI = { 2,   /* nDims */
  115,                                 /* lineNo */
  14,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtDCInfo d_emlrtDCI = { 125, /* lineNo */
  79,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  1                                    /* checkKind */
};

static emlrtECInfo d_emlrtECI = { 2,   /* nDims */
  125,                                 /* lineNo */
  50,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtDCInfo e_emlrtDCI = { 127, /* lineNo */
  30,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  1                                    /* checkKind */
};

static emlrtECInfo e_emlrtECI = { 2,   /* nDims */
  127,                                 /* lineNo */
  16,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtECInfo f_emlrtECI = { 2,   /* nDims */
  128,                                 /* lineNo */
  29,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtBCInfo g_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  97,                                  /* lineNo */
  26,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  142,                                 /* lineNo */
  22,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtECInfo g_emlrtECI = { -1,  /* nDims */
  142,                                 /* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtBCInfo i_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  142,                                 /* lineNo */
  39,                                  /* colNo */
  "mt_spectrum",                       /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo b_emlrtRTEI = { 20,/* lineNo */
  15,                                  /* colNo */
  "rdivide_helper",                    /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/rdivide_helper.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 283,/* lineNo */
  27,                                  /* colNo */
  "check_non_axis_size",               /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/cat.m"/* pName */
};

static emlrtBCInfo j_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  121,                                 /* lineNo */
  37,                                  /* colNo */
  "Spower",                            /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m",         /* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo w_emlrtRTEI = { 55,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo x_emlrtRTEI = { 65,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo y_emlrtRTEI = { 68,/* lineNo */
  14,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = { 68,/* lineNo */
  47,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = { 68,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = { 1,/* lineNo */
  45,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo db_emlrtRTEI = { 73,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = { 68,/* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = { 77,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = { 146,/* lineNo */
  18,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = { 146,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = { 147,/* lineNo */
  23,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo jb_emlrtRTEI = { 28,/* lineNo */
  9,                                   /* colNo */
  "colon",                             /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/colon.m"/* pName */
};

static emlrtRTEInfo kb_emlrtRTEI = { 148,/* lineNo */
  36,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = { 89,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo mb_emlrtRTEI = { 149,/* lineNo */
  50,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo nb_emlrtRTEI = { 92,/* lineNo */
  12,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = { 149,/* lineNo */
  19,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = { 17,/* lineNo */
  13,                                  /* colNo */
  "isnan",                             /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elmat/isnan.m"/* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = { 149,/* lineNo */
  76,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = { 103,/* lineNo */
  32,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo sb_emlrtRTEI = { 149,/* lineNo */
  34,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo tb_emlrtRTEI = { 105,/* lineNo */
  32,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo ub_emlrtRTEI = { 149,/* lineNo */
  91,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = { 115,/* lineNo */
  14,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = { 115,/* lineNo */
  34,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = { 137,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo yb_emlrtRTEI = { 133,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo ac_emlrtRTEI = { 138,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo bc_emlrtRTEI = { 134,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo cc_emlrtRTEI = { 121,/* lineNo */
  28,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo dc_emlrtRTEI = { 149,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo ec_emlrtRTEI = { 122,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo fc_emlrtRTEI = { 142,/* lineNo */
  27,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo gc_emlrtRTEI = { 125,/* lineNo */
  50,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo hc_emlrtRTEI = { 125,/* lineNo */
  74,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo ic_emlrtRTEI = { 99,/* lineNo */
  9,                                   /* colNo */
  "eml_mtimes_helper",                 /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pName */
};

static emlrtRTEInfo jc_emlrtRTEI = { 127,/* lineNo */
  25,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo kc_emlrtRTEI = { 128,/* lineNo */
  29,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo lc_emlrtRTEI = { 128,/* lineNo */
  34,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo mc_emlrtRTEI = { 128,/* lineNo */
  25,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo nc_emlrtRTEI = { 127,/* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo oc_emlrtRTEI = { 125,/* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo pc_emlrtRTEI = { 121,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo qc_emlrtRTEI = { 115,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo rc_emlrtRTEI = { 112,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo sc_emlrtRTEI = { 89,/* lineNo */
  44,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo tc_emlrtRTEI = { 128,/* lineNo */
  47,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo uc_emlrtRTEI = { 142,/* lineNo */
  39,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo vc_emlrtRTEI = { 148,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

static emlrtRTEInfo wc_emlrtRTEI = { 33,/* lineNo */
  6,                                   /* colNo */
  "find",                              /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elmat/find.m"/* pName */
};

static emlrtRTEInfo xc_emlrtRTEI = { 148,/* lineNo */
  18,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/Prerau_labcode/multitaper/multitaper_spectrogram_mex/multi"
  "taper_spectrogram_coder.m"          /* pName */
};

/* Function Declarations */
static void b_dynamic_size_checks(const emlrtStack *sp, const emxArray_real32_T *
  a, const emxArray_real_T *b, int32_T innerDimA, int32_T innerDimB);
static void dynamic_size_checks(const emlrtStack *sp, const emxArray_real32_T *a,
  const emxArray_real32_T *b, int32_T innerDimA, int32_T innerDimB);

/* Function Definitions */
static void b_dynamic_size_checks(const emlrtStack *sp, const emxArray_real32_T *
  a, const emxArray_real_T *b, int32_T innerDimA, int32_T innerDimB)
{
  if (innerDimA != innerDimB) {
    if (((a->size[0] == 1) && (a->size[1] == 1)) || (b->size[0] == 1)) {
      emlrtErrorWithMessageIdR2018a(sp, &j_emlrtRTEI,
        "Coder:toolbox:mtimes_noDynamicScalarExpansion",
        "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
    } else {
      emlrtErrorWithMessageIdR2018a(sp, &k_emlrtRTEI, "MATLAB:innerdim",
        "MATLAB:innerdim", 0);
    }
  }
}

static void dynamic_size_checks(const emlrtStack *sp, const emxArray_real32_T *a,
  const emxArray_real32_T *b, int32_T innerDimA, int32_T innerDimB)
{
  if (innerDimA != innerDimB) {
    if ((a->size[0] == 1) || (b->size[0] == 1)) {
      emlrtErrorWithMessageIdR2018a(sp, &j_emlrtRTEI,
        "Coder:toolbox:mtimes_noDynamicScalarExpansion",
        "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
    } else {
      emlrtErrorWithMessageIdR2018a(sp, &k_emlrtRTEI, "MATLAB:innerdim",
        "MATLAB:innerdim", 0);
    }
  }
}

void multitaper_spectrogram_coder(const emlrtStack *sp, const emxArray_real32_T *
  data, real_T Fs, const real_T frequency_range[2], const emxArray_real_T
  *DPSS_tapers, const emxArray_real_T *DPSS_eigen, real_T winstep_samples,
  real_T min_NFFT, real_T detrend_opt, real_T weighting, emxArray_real32_T
  *mt_spectrogram, emxArray_real_T *stimes, emxArray_real_T *sfreqs)
{
  int32_T num_tapers;
  int32_T winsize_samples;
  real_T b;
  emxArray_real_T *window_start;
  int32_T i;
  int32_T loop_ub;
  real_T nfft;
  emxArray_boolean_T *x;
  emxArray_boolean_T *r;
  emxArray_boolean_T *freq_inds;
  int32_T nz;
  int32_T trueCount;
  int32_T b_i;
  int32_T input_sizes_idx_0;
  int32_T n;
  emxArray_int32_T *r1;
  emxArray_int32_T *r2;
  emxArray_int32_T *r3;
  emxArray_real_T *r4;
  emxArray_real32_T *y;
  emxArray_real32_T *b_y;
  emxArray_real32_T *Spower;
  emxArray_real32_T *data_segment;
  emxArray_boolean_T *r5;
  emxArray_boolean_T *b_data_segment;
  emxArray_int32_T *ii;
  emxArray_real_T *c_y;
  emxArray_real32_T *c_data_segment;
  emxArray_creal32_T *fft_data;
  emxArray_real32_T *b_Spower;
  emxArray_real32_T *a;
  emxArray_real32_T *Spower_iter;
  emxArray_real32_T *b_b;
  emxArray_real32_T *wk;
  emxArray_real_T *wt;
  emxArray_real_T *d_y;
  real32_T Tpower;
  emxArray_real_T *b_window_start;
  emxArray_real_T *b_select;
  emxArray_real32_T *c_b;
  emxArray_int32_T *ia;
  int32_T ib_size[1];
  int32_T i1;
  int32_T b_loop_ub;
  real_T c_window_start;
  real_T d;
  emxArray_real32_T *varargin_1;
  int32_T end;
  emxArray_real32_T *varargin_3;
  int32_T d_window_start[1];
  emxArray_real32_T d_data_segment;
  int32_T e_window_start[1];
  boolean_T empty_non_axis_sizes;
  boolean_T guard1 = false;
  int32_T f_window_start[1];
  emxArray_real32_T e_data_segment;
  int32_T f_data_segment[1];
  int32_T c_loop_ub;
  int32_T g_window_start[1];
  real32_T f;
  int32_T g_data_segment[1];
  int32_T d_loop_ub;
  int32_T e_loop_ub;
  int32_T b_ii;
  int32_T f_loop_ub;
  int32_T i2;
  int32_T g_loop_ub;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  boolean_T emlrtHadParallelError = false;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /* MULTITAPER_SPECTROGRAM  Compute the multitaper spectrogram for time series data */
  /*  */
  /*    This is the coder portion for mex compilation. It takes processed */
  /*    multitaper parameters from multitaper_spectrogram_coder_mex.m as inputs */
  /*    and skips internal input processing.  */
  /*     */
  /*    Usage: */
  /*    Direct input: */
  /*    [spect,stimes,sfreqs] = multitaper_spectrogram_coder(data, Fs, frequency_range, DPSS_tapers, DPSS_eigen, winstep_samples, min_NFFT, detrend_opt, weighting) */
  /*  */
  /*    Input: */
  /*    data: 1 x <number of samples> vector - time series data-- required */
  /*    Fs: double - sampling frequency in Hz  -- required */
  /*    frequency_range: 1x2 vector - [<min frequency>, <max frequency>] -- required  */
  /*    DPSS_tapers: Nxk matrix - Slepian tapers -- required  */
  /*    DPSS_eigen: 1xk vector - eigenvalues of the Slepian tapers -- required */
  /*    winstep_samples: double - number of samples in step size of windows -- required  */
  /*    min_NFFT: double - minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) -- required */
  /*    detrend_opt: double - how to detrend data window (2='linear', 1='constant', 0='off') -- required */
  /*    weighting: double - how to weight the tapers (0='unity', 1='eigen', 2='adapt') -- required */
  /*  */
  /*    Output: */
  /*    mt_spectrogram: TxF matrix of one-sided power spectral density (PSD) */
  /*    stimes: 1XT vector of times for the center of the spectral bins */
  /*    sfreqs: 1XF vector of frequency bins for the spectrogram */
  /*  */
  /*    This code is companion to the paper: */
  /*          "Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis" */
  /*          Michael J. Prerau, Ritchie E. Brown, Matt T. Bianchi, Jeffrey M. Ellenbogen, Patrick L. Purdon */
  /*          December 7, 2016 : 60-92 */
  /*          DOI: 10.1152/physiol.00062.2015 */
  /*    which should be cited for academic use of this code. */
  /*  */
  /*    A full tutorial on the multitaper spectrogram can be found at: */
  /*    http://www.sleepEEG.org/multitaper */
  /*  */
  /*    Copyright 2019 Michael J. Prerau, Ph.D. - http://www.sleepEEG.org */
  /*    This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. */
  /*    (http://creativecommons.org/licenses/by-nc-sa/4.0/) */
  /*  */
  /*    Last modified 2/16/2021 */
  /*  ******************************************************************** */
  /* Generate DPSS tapers (STEP 1) */
  /*  Done outside this _coder function. */
  /* Get taper matrix dimensions  */
  num_tapers = DPSS_tapers->size[1];
  winsize_samples = DPSS_tapers->size[0];

  /* Total data length */
  /* Window start indices */
  st.site = &emlrtRSI;
  b = static_cast<real_T>(data->size[1] - DPSS_tapers->size[0]) + 1.0;
  emxInit_real_T(&st, &window_start, 2, &w_emlrtRTEI, true);
  if (muDoubleScalarIsNaN(winstep_samples)) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &w_emlrtRTEI);
    window_start->data[0] = rtNaN;
  } else if ((winstep_samples == 0.0) || ((1.0 < b) && (winstep_samples < 0.0)) ||
             ((b < 1.0) && (winstep_samples > 0.0))) {
    window_start->size[0] = 1;
    window_start->size[1] = 0;
  } else if (muDoubleScalarIsInf(winstep_samples)) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &w_emlrtRTEI);
    window_start->data[0] = 1.0;
  } else if (muDoubleScalarFloor(winstep_samples) == winstep_samples) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    loop_ub = static_cast<int32_T>(muDoubleScalarFloor((b - 1.0) /
      winstep_samples));
    window_start->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &w_emlrtRTEI);
    for (i = 0; i <= loop_ub; i++) {
      window_start->data[i] = winstep_samples * static_cast<real_T>(i) + 1.0;
    }
  } else {
    b_st.site = &w_emlrtRSI;
    eml_float_colon(&b_st, winstep_samples, b, window_start);
  }

  /* Number of windows */
  /* Number of points in the FFT */
  st.site = &b_emlrtRSI;
  b_st.site = &b_emlrtRSI;
  b = nextpow2(static_cast<real_T>(DPSS_tapers->size[0]));
  b_st.site = &gb_emlrtRSI;
  st.site = &b_emlrtRSI;
  b_st.site = &b_emlrtRSI;
  nfft = nextpow2(min_NFFT);
  b_st.site = &gb_emlrtRSI;
  nfft = muDoubleScalarMax(muDoubleScalarMax(muDoubleScalarPower(2.0, b),
    static_cast<real_T>(DPSS_tapers->size[0])), muDoubleScalarPower(2.0, nfft));

  /* Create the frequency vector */
  b = Fs / nfft;
  st.site = &c_emlrtRSI;
  if (muDoubleScalarIsNaN(b) || muDoubleScalarIsNaN(Fs)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &x_emlrtRTEI);
    sfreqs->data[0] = rtNaN;
  } else if ((b == 0.0) || ((0.0 < Fs) && (b < 0.0)) || ((Fs < 0.0) && (b > 0.0)))
  {
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 0;
  } else if (muDoubleScalarIsInf(Fs) && (muDoubleScalarIsInf(b) || (0.0 == Fs)))
  {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &x_emlrtRTEI);
    sfreqs->data[0] = rtNaN;
  } else if (muDoubleScalarIsInf(b)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &x_emlrtRTEI);
    sfreqs->data[0] = 0.0;
  } else if (muDoubleScalarFloor(b) == b) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    loop_ub = static_cast<int32_T>(muDoubleScalarFloor(Fs / b));
    sfreqs->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &x_emlrtRTEI);
    for (i = 0; i <= loop_ub; i++) {
      sfreqs->data[i] = b * static_cast<real_T>(i);
    }
  } else {
    b_st.site = &w_emlrtRSI;
    b_eml_float_colon(&b_st, b, Fs, sfreqs);
  }

  emxInit_boolean_T(&st, &x, 2, &eb_emlrtRTEI, true);

  /*  all possible frequencies */
  /* Get just the frequencies for the given frequency range */
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, x, i, &y_emlrtRTEI);
  loop_ub = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < loop_ub; i++) {
    x->data[i] = (sfreqs->data[i] >= frequency_range[0]);
  }

  emxInit_boolean_T(sp, &r, 2, &cb_emlrtRTEI, true);
  i = r->size[0] * r->size[1];
  r->size[0] = 1;
  r->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, r, i, &ab_emlrtRTEI);
  loop_ub = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < loop_ub; i++) {
    r->data[i] = (sfreqs->data[i] <= frequency_range[1]);
  }

  emxInit_boolean_T(sp, &freq_inds, 2, &bb_emlrtRTEI, true);
  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])x->size, *(int32_T (*)[2])r->size,
    &emlrtECI, sp);
  i = freq_inds->size[0] * freq_inds->size[1];
  freq_inds->size[0] = 1;
  freq_inds->size[1] = x->size[1];
  emxEnsureCapacity_boolean_T(sp, freq_inds, i, &bb_emlrtRTEI);
  loop_ub = x->size[0] * x->size[1];
  for (i = 0; i < loop_ub; i++) {
    freq_inds->data[i] = (x->data[i] && r->data[i]);
  }

  nz = x->size[1] - 1;
  trueCount = 0;
  for (b_i = 0; b_i <= nz; b_i++) {
    if (x->data[b_i] && r->data[b_i]) {
      trueCount++;
    }
  }

  input_sizes_idx_0 = 0;
  for (b_i = 0; b_i <= nz; b_i++) {
    if (x->data[b_i] && r->data[b_i]) {
      i = b_i + 1;
      if ((i < 1) || (i > sfreqs->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i, 1, sfreqs->size[1], &f_emlrtBCI, sp);
      }

      sfreqs->data[input_sizes_idx_0] = sfreqs->data[i - 1];
      input_sizes_idx_0++;
    }
  }

  i = sfreqs->size[0] * sfreqs->size[1];
  sfreqs->size[0] = 1;
  sfreqs->size[1] = trueCount;
  emxEnsureCapacity_real_T(sp, sfreqs, i, &cb_emlrtRTEI);

  /* Compute the times of the middle of each spectrum */
  nz = static_cast<int32_T>(muDoubleScalarRound(static_cast<real_T>
    (DPSS_tapers->size[0]) / 2.0));
  i = stimes->size[0] * stimes->size[1];
  stimes->size[0] = 1;
  stimes->size[1] = window_start->size[1];
  emxEnsureCapacity_real_T(sp, stimes, i, &db_emlrtRTEI);
  loop_ub = window_start->size[0] * window_start->size[1];
  for (i = 0; i < loop_ub; i++) {
    stimes->data[i] = (window_start->data[i] + static_cast<real_T>(nz)) / Fs;
  }

  /* Preallocate spectrogram and slice data for efficient parallel computing */
  st.site = &d_emlrtRSI;
  i = x->size[0] * x->size[1];
  input_sizes_idx_0 = x->size[0] * x->size[1];
  x->size[0] = 1;
  emxEnsureCapacity_boolean_T(&st, x, input_sizes_idx_0, &eb_emlrtRTEI);
  loop_ub = i - 1;
  for (i = 0; i <= loop_ub; i++) {
    x->data[i] = (x->data[i] && r->data[i]);
  }

  emxFree_boolean_T(&r);
  b_st.site = &ib_emlrtRSI;
  c_st.site = &jb_emlrtRSI;
  nz = combineVectorElements(&c_st, x);
  if (nz < 0) {
    emlrtNonNegativeCheckR2012b(static_cast<real_T>(nz), &c_emlrtDCI, sp);
  }

  i = mt_spectrogram->size[0] * mt_spectrogram->size[1];
  mt_spectrogram->size[0] = nz;
  mt_spectrogram->size[1] = window_start->size[1];
  emxEnsureCapacity_real32_T(sp, mt_spectrogram, i, &fb_emlrtRTEI);
  loop_ub = nz * window_start->size[1];
  for (i = 0; i < loop_ub; i++) {
    mt_spectrogram->data[i] = 0.0F;
  }

  nz = window_start->size[1] - 1;
  emlrtEnterParallelRegion(sp, omp_in_parallel());

#pragma omp parallel \
 num_threads(emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs())) \
 private(r1,r2,r3,r4,y,b_y,Spower,data_segment,r5,b_data_segment,c_y,c_data_segment,fft_data,b_Spower,a,Spower_iter,b_b,wk,wt,Tpower,i_st,i1,c_window_start,b_loop_ub,d,end,d_window_start,d_data_segment,e_data_segment,f,c_loop_ub,d_loop_ub,e_loop_ub,b_ii,f_loop_ub,i2,g_loop_ub) \
 firstprivate(d_st,e_st,f_st,g_st,h_st,e_window_start,f_window_start,f_data_segment,emlrtHadParallelError,g_window_start,g_data_segment)

  {
    try {
      d_st.prev = sp;
      d_st.tls = emlrtAllocTLS((emlrtStack *)sp, omp_get_thread_num());
      d_st.site = NULL;
      e_st.prev = &d_st;
      e_st.tls = d_st.tls;
      f_st.prev = &e_st;
      f_st.tls = e_st.tls;
      g_st.prev = &f_st;
      g_st.tls = f_st.tls;
      h_st.prev = &g_st;
      h_st.tls = g_st.tls;
      i_st.prev = &h_st;
      i_st.tls = h_st.tls;
      emxInit_int32_T(&d_st, &r1, 2, &uc_emlrtRTEI, true);
      emxInit_int32_T(&d_st, &r2, 1, &fc_emlrtRTEI, true);
      emxInit_int32_T(&d_st, &r3, 2, &uc_emlrtRTEI, true);
      emxInit_real_T(&d_st, &r4, 2, &cb_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &y, 2, &gc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &b_y, 1, &tc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &Spower, 2, &cc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &data_segment, 2, &tb_emlrtRTEI, true);
      emxInit_boolean_T(&d_st, &r5, 2, &pb_emlrtRTEI, true);
      emxInit_boolean_T(&d_st, &b_data_segment, 2, &nb_emlrtRTEI, true);
      emxInit_real_T(&d_st, &c_y, 2, &sc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &c_data_segment, 2, &lb_emlrtRTEI, true);
      emxInit_creal32_T(&d_st, &fft_data, 2, &rc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &b_Spower, 2, &qc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &a, 1, &ec_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &Spower_iter, 1, &pc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &b_b, 2, &oc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &wk, 2, &nc_emlrtRTEI, true);
      emxInit_real_T(&d_st, &wt, 1, &yb_emlrtRTEI, true);
    } catch (...) {
      emlrtHadParallelError = true;
    }

#pragma omp for nowait

    for (n = 0; n <= nz; n++) {
      if (emlrtHadParallelError)
        continue;
      try {
        /*  COMPUTE THE MULTITAPER SPECTROGRAM */
        /*  */
        /*      STEP 1: Compute DPSS tapers based on desired spectral properties */
        /*      STEP 2: Multiply the data segment by the DPSS Tapers */
        /*      STEP 3: Compute the spectrum for each tapered segment */
        /*      STEP 4: Take the mean of the tapered spectra */
        /* Loop in parallel over all of the windows */
        /* Grab the data for the given window */
        if (winsize_samples - 1 < 0) {
          c_y->size[0] = 1;
          c_y->size[1] = 0;
        } else {
          i1 = c_y->size[0] * c_y->size[1];
          c_y->size[0] = 1;
          b_loop_ub = static_cast<int32_T>(static_cast<real_T>(winsize_samples)
            - 1.0);
          c_y->size[1] = b_loop_ub + 1;
          emxEnsureCapacity_real_T(&d_st, c_y, i1, &jb_emlrtRTEI);
          for (i1 = 0; i1 <= b_loop_ub; i1++) {
            c_y->data[i1] = i1;
          }
        }

        i1 = c_data_segment->size[0] * c_data_segment->size[1];
        c_data_segment->size[0] = 1;
        c_data_segment->size[1] = c_y->size[1];
        emxEnsureCapacity_real32_T(&d_st, c_data_segment, i1, &lb_emlrtRTEI);
        i1 = n + 1;
        if ((i1 < 1) || (i1 > window_start->size[1])) {
          emlrtDynamicBoundsCheckR2012b(i1, 1, window_start->size[1],
            &d_emlrtBCI, &d_st);
        }

        c_window_start = window_start->data[i1 - 1];
        b_loop_ub = c_y->size[0] * c_y->size[1];
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          d = c_window_start + c_y->data[i1];
          if (d != static_cast<int32_T>(muDoubleScalarFloor(d))) {
            emlrtIntegerCheckR2012b(d, &emlrtDCI, &d_st);
          }

          end = static_cast<int32_T>(d);
          if ((end < 1) || (end > data->size[1])) {
            emlrtDynamicBoundsCheckR2012b(end, 1, data->size[1], &c_emlrtBCI,
              &d_st);
          }

          c_data_segment->data[i1] = data->data[end - 1];
        }

        /* Skip empty segments */
        i1 = b_data_segment->size[0] * b_data_segment->size[1];
        b_data_segment->size[0] = 1;
        b_data_segment->size[1] = c_data_segment->size[1];
        emxEnsureCapacity_boolean_T(&d_st, b_data_segment, i1, &nb_emlrtRTEI);
        b_loop_ub = c_data_segment->size[0] * c_data_segment->size[1];
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          b_data_segment->data[i1] = (c_data_segment->data[i1] == 0.0F);
        }

        e_st.site = &e_emlrtRSI;
        if (!all(&e_st, b_data_segment)) {
          i1 = r5->size[0] * r5->size[1];
          r5->size[0] = 1;
          r5->size[1] = c_data_segment->size[1];
          emxEnsureCapacity_boolean_T(&d_st, r5, i1, &pb_emlrtRTEI);
          b_loop_ub = c_data_segment->size[0] * c_data_segment->size[1];
          for (i1 = 0; i1 < b_loop_ub; i1++) {
            r5->data[i1] = muSingleScalarIsNaN(c_data_segment->data[i1]);
          }

          e_st.site = &f_emlrtRSI;
          if (any(&e_st, r5)) {
            i1 = n + 1;
            if ((i1 < 1) || (i1 > mt_spectrogram->size[1])) {
              emlrtDynamicBoundsCheckR2012b(i1, 1, mt_spectrogram->size[1],
                &g_emlrtBCI, &d_st);
            }

            b_loop_ub = mt_spectrogram->size[0];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              mt_spectrogram->data[i1 + mt_spectrogram->size[0] * n] = rtNaNF;
            }
          } else {
            /* Option to detrend_opt data to remove low frequency DC component */
            if (detrend_opt == 1.0) {
              i1 = data_segment->size[0] * data_segment->size[1];
              data_segment->size[0] = 1;
              data_segment->size[1] = c_data_segment->size[1];
              emxEnsureCapacity_real32_T(&d_st, data_segment, i1, &rb_emlrtRTEI);
              b_loop_ub = c_data_segment->size[0] * c_data_segment->size[1] - 1;
              for (i1 = 0; i1 <= b_loop_ub; i1++) {
                data_segment->data[i1] = c_data_segment->data[i1];
              }

              e_st.site = &g_emlrtRSI;
              detrend(&e_st, data_segment, c_data_segment);
            } else {
              if (detrend_opt == 2.0) {
                i1 = data_segment->size[0] * data_segment->size[1];
                data_segment->size[0] = 1;
                data_segment->size[1] = c_data_segment->size[1];
                emxEnsureCapacity_real32_T(&d_st, data_segment, i1,
                  &tb_emlrtRTEI);
                b_loop_ub = c_data_segment->size[0] * c_data_segment->size[1] -
                  1;
                for (i1 = 0; i1 <= b_loop_ub; i1++) {
                  data_segment->data[i1] = c_data_segment->data[i1];
                }

                e_st.site = &h_emlrtRSI;
                b_detrend(&e_st, data_segment, c_data_segment);
              }
            }

            /* Multiply the data by the tapers (STEP 2) */
            d_window_start[0] = c_data_segment->size[1];
            d_data_segment = *c_data_segment;
            e_window_start[0] = d_window_start[0];
            d_data_segment.size = &e_window_start[0];
            d_data_segment.numDimensions = 1;
            e_st.site = &i_emlrtRSI;
            repmat(&e_st, &d_data_segment, static_cast<real_T>(num_tapers), wk);
            emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])wk->size, *(int32_T (*)[2])
              DPSS_tapers->size, &b_emlrtECI, &d_st);
            b_loop_ub = wk->size[0] * wk->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              wk->data[i1] *= static_cast<real32_T>(DPSS_tapers->data[i1]);
            }

            /* Compute the FFT (STEP 3) */
            e_st.site = &j_emlrtRSI;
            fft(&e_st, wk, nfft, fft_data);

            /* Compute the weighted mean spectral power across tapers (STEP 4) */
            i1 = b_b->size[0] * b_b->size[1];
            b_b->size[0] = fft_data->size[0];
            b_b->size[1] = fft_data->size[1];
            emxEnsureCapacity_real32_T(&d_st, b_b, i1, &vb_emlrtRTEI);
            b_loop_ub = fft_data->size[0] * fft_data->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              b_b->data[i1] = fft_data->data[i1].im;
            }

            e_st.site = &k_emlrtRSI;
            power(&e_st, b_b, b_Spower);
            i1 = b_b->size[0] * b_b->size[1];
            b_b->size[0] = fft_data->size[0];
            b_b->size[1] = fft_data->size[1];
            emxEnsureCapacity_real32_T(&d_st, b_b, i1, &wb_emlrtRTEI);
            b_loop_ub = fft_data->size[0] * fft_data->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              b_b->data[i1] = fft_data->data[i1].re;
            }

            e_st.site = &k_emlrtRSI;
            power(&e_st, b_b, wk);
            emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])b_Spower->size, *(int32_T
              (*)[2])wk->size, &c_emlrtECI, &d_st);
            b_loop_ub = b_Spower->size[0] * b_Spower->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              b_Spower->data[i1] += wk->data[i1];
            }

            if (weighting == 2.0) {
              /*  adaptive weights - for colored noise spectrum (Percival & Walden */
              /*  p368-p370) */
              e_st.site = &l_emlrtRSI;
              d_window_start[0] = c_data_segment->size[1];
              b_loop_ub = c_data_segment->size[1];
              d_data_segment = *c_data_segment;
              f_window_start[0] = d_window_start[0];
              d_data_segment.size = &f_window_start[0];
              d_data_segment.numDimensions = 1;
              e_data_segment = *c_data_segment;
              f_data_segment[0] = b_loop_ub;
              e_data_segment.size = &f_data_segment[0];
              e_data_segment.numDimensions = 1;
              f_st.site = &fe_emlrtRSI;
              dynamic_size_checks(&f_st, &d_data_segment, &e_data_segment,
                                  c_data_segment->size[1], c_data_segment->size
                                  [1]);
              if (c_data_segment->size[1] == 1) {
                Tpower = 0.0F;
                b_loop_ub = c_data_segment->size[1];
                for (i1 = 0; i1 < b_loop_ub; i1++) {
                  f = c_data_segment->data[i1];
                  Tpower += f * f;
                }
              } else {
                d_window_start[0] = c_data_segment->size[1];
                b_loop_ub = c_data_segment->size[1];
                d_data_segment = *c_data_segment;
                g_window_start[0] = d_window_start[0];
                d_data_segment.size = &g_window_start[0];
                d_data_segment.numDimensions = 1;
                e_data_segment = *c_data_segment;
                g_data_segment[0] = b_loop_ub;
                e_data_segment.size = &g_data_segment[0];
                e_data_segment.numDimensions = 1;
                f_st.site = &ee_emlrtRSI;
                Tpower = mtimes(&f_st, &d_data_segment, &e_data_segment);
              }

              Tpower /= static_cast<real32_T>(c_data_segment->size[1]);
              b_loop_ub = b_Spower->size[0];
              i1 = Spower->size[0] * Spower->size[1];
              Spower->size[0] = b_Spower->size[0];
              Spower->size[1] = 2;
              emxEnsureCapacity_real32_T(&d_st, Spower, i1, &cc_emlrtRTEI);
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                if (1 > b_Spower->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(1, 1, b_Spower->size[1],
                    &j_emlrtBCI, &d_st);
                }

                Spower->data[i1] = b_Spower->data[i1];
              }

              for (i1 = 0; i1 < b_loop_ub; i1++) {
                if (2 > b_Spower->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(2, 1, b_Spower->size[1],
                    &j_emlrtBCI, &d_st);
                }

                Spower->data[i1 + Spower->size[0]] = b_Spower->data[i1 +
                  b_Spower->size[0]];
              }

              e_st.site = &m_emlrtRSI;
              mean(&e_st, Spower, Spower_iter);
              b_loop_ub = DPSS_eigen->size[0];
              i1 = a->size[0];
              a->size[0] = DPSS_eigen->size[0];
              emxEnsureCapacity_real32_T(&d_st, a, i1, &ec_emlrtRTEI);
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                a->data[i1] = static_cast<real32_T>(1.0 - DPSS_eigen->data[i1]) *
                  Tpower;
              }

              b_loop_ub = DPSS_eigen->size[0];
              c_loop_ub = a->size[0];
              d_loop_ub = b_Spower->size[0];
              i1 = static_cast<int32_T>(muDoubleScalarFloor(nfft));
              e_loop_ub = static_cast<int32_T>(nfft);
              for (b_ii = 0; b_ii < 3; b_ii++) {
                /*  run 3 iterations */
                /*  calculate the MSE weights */
                end = y->size[0] * y->size[1];
                y->size[0] = Spower_iter->size[0];
                y->size[1] = DPSS_eigen->size[0];
                emxEnsureCapacity_real32_T(&d_st, y, end, &gc_emlrtRTEI);
                for (end = 0; end < b_loop_ub; end++) {
                  f_loop_ub = Spower_iter->size[0];
                  for (i2 = 0; i2 < f_loop_ub; i2++) {
                    y->data[i2 + y->size[0] * end] = Spower_iter->data[i2] *
                      static_cast<real32_T>(DPSS_eigen->data[end]);
                  }
                }

                if (nfft != i1) {
                  emlrtIntegerCheckR2012b(nfft, &d_emlrtDCI, &d_st);
                }

                end = wk->size[0] * wk->size[1];
                wk->size[0] = static_cast<int32_T>(nfft);
                wk->size[1] = a->size[0];
                emxEnsureCapacity_real32_T(&d_st, wk, end, &hc_emlrtRTEI);
                for (end = 0; end < c_loop_ub; end++) {
                  for (i2 = 0; i2 < e_loop_ub; i2++) {
                    wk->data[i2 + wk->size[0] * end] = a->data[end];
                  }
                }

                emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])y->size, *(int32_T (*)
                  [2])wk->size, &d_emlrtECI, &d_st);
                end = b_b->size[0] * b_b->size[1];
                b_b->size[0] = Spower_iter->size[0];
                b_b->size[1] = num_tapers;
                emxEnsureCapacity_real32_T(&d_st, b_b, end, &ic_emlrtRTEI);
                for (end = 0; end < num_tapers; end++) {
                  f_loop_ub = Spower_iter->size[0];
                  for (i2 = 0; i2 < f_loop_ub; i2++) {
                    b_b->data[i2 + b_b->size[0] * end] = Spower_iter->data[i2];
                  }
                }

                e_st.site = &n_emlrtRSI;
                f_loop_ub = y->size[0] * y->size[1];
                for (end = 0; end < f_loop_ub; end++) {
                  y->data[end] += wk->data[end];
                }

                if (!scalexpCompatible(b_b, y)) {
                  emlrtErrorWithMessageIdR2018a(&e_st, &b_emlrtRTEI,
                    "MATLAB:dimagree", "MATLAB:dimagree", 0);
                }

                f_loop_ub = b_b->size[0] * b_b->size[1];
                for (end = 0; end < f_loop_ub; end++) {
                  b_b->data[end] /= y->data[end];
                }

                /*  calculate new spectral estimate */
                e_st.site = &o_emlrtRSI;
                power(&e_st, b_b, wk);
                if (static_cast<int32_T>(nfft) != i1) {
                  emlrtIntegerCheckR2012b(nfft, &e_emlrtDCI, &d_st);
                }

                end = r4->size[0] * r4->size[1];
                r4->size[0] = static_cast<int32_T>(nfft);
                r4->size[1] = DPSS_eigen->size[0];
                emxEnsureCapacity_real_T(&d_st, r4, end, &jc_emlrtRTEI);
                for (end = 0; end < b_loop_ub; end++) {
                  for (i2 = 0; i2 < e_loop_ub; i2++) {
                    r4->data[i2 + r4->size[0] * end] = DPSS_eigen->data[end];
                  }
                }

                emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])wk->size, *(int32_T (*)
                  [2])r4->size, &e_emlrtECI, &d_st);
                f_loop_ub = wk->size[0] * wk->size[1];
                for (end = 0; end < f_loop_ub; end++) {
                  wk->data[end] *= static_cast<real32_T>(r4->data[end]);
                }

                end = b_b->size[0] * b_b->size[1];
                b_b->size[0] = wk->size[1];
                f_loop_ub = wk->size[0];
                b_b->size[1] = wk->size[0];
                emxEnsureCapacity_real32_T(&d_st, b_b, end, &kc_emlrtRTEI);
                for (end = 0; end < f_loop_ub; end++) {
                  g_loop_ub = wk->size[1];
                  for (i2 = 0; i2 < g_loop_ub; i2++) {
                    b_b->data[i2 + b_b->size[0] * end] = wk->data[end + wk->
                      size[0] * i2];
                  }
                }

                end = y->size[0] * y->size[1];
                y->size[0] = b_Spower->size[1];
                y->size[1] = b_Spower->size[0];
                emxEnsureCapacity_real32_T(&d_st, y, end, &lc_emlrtRTEI);
                for (end = 0; end < d_loop_ub; end++) {
                  f_loop_ub = b_Spower->size[1];
                  for (i2 = 0; i2 < f_loop_ub; i2++) {
                    y->data[i2 + y->size[0] * end] = b_Spower->data[end +
                      b_Spower->size[0] * i2];
                  }
                }

                emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])b_b->size, *(int32_T (*)
                  [2])y->size, &f_emlrtECI, &d_st);
                e_st.site = &p_emlrtRSI;
                end = b_b->size[0] * b_b->size[1];
                b_b->size[0] = wk->size[1];
                f_loop_ub = wk->size[0];
                b_b->size[1] = wk->size[0];
                emxEnsureCapacity_real32_T(&e_st, b_b, end, &kc_emlrtRTEI);
                for (end = 0; end < f_loop_ub; end++) {
                  g_loop_ub = wk->size[1];
                  for (i2 = 0; i2 < g_loop_ub; i2++) {
                    b_b->data[i2 + b_b->size[0] * end] = wk->data[end + wk->
                      size[0] * i2] * b_Spower->data[end + b_Spower->size[0] *
                      i2];
                  }
                }

                f_st.site = &p_emlrtRSI;
                sum(&f_st, b_b, c_data_segment);
                end = Spower_iter->size[0];
                Spower_iter->size[0] = c_data_segment->size[1];
                emxEnsureCapacity_real32_T(&e_st, Spower_iter, end,
                  &mc_emlrtRTEI);
                f_loop_ub = c_data_segment->size[1];
                for (end = 0; end < f_loop_ub; end++) {
                  Spower_iter->data[end] = c_data_segment->data[end];
                }

                f_st.site = &p_emlrtRSI;
                g_st.site = &lf_emlrtRSI;
                h_st.site = &jb_emlrtRSI;
                i_st.site = &hf_emlrtRSI;
                blockedSummation(&i_st, wk, wk->size[1], b_y);
                if (!b_scalexpCompatible(Spower_iter, b_y)) {
                  emlrtErrorWithMessageIdR2018a(&e_st, &b_emlrtRTEI,
                    "MATLAB:dimagree", "MATLAB:dimagree", 0);
                }

                f_loop_ub = Spower_iter->size[0];
                for (end = 0; end < f_loop_ub; end++) {
                  Spower_iter->data[end] /= b_y->data[end];
                }

                if (*emlrtBreakCheckR2012bFlagVar != 0) {
                  emlrtBreakCheckR2012b(&d_st);
                }
              }
            } else if (weighting == 1.0) {
              /*  eigenvalue weights */
              b_loop_ub = DPSS_eigen->size[0];
              i1 = wt->size[0];
              wt->size[0] = DPSS_eigen->size[0];
              emxEnsureCapacity_real_T(&d_st, wt, i1, &yb_emlrtRTEI);
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                wt->data[i1] = DPSS_eigen->data[i1] / static_cast<real_T>
                  (num_tapers);
              }

              e_st.site = &q_emlrtRSI;
              f_st.site = &fe_emlrtRSI;
              b_dynamic_size_checks(&f_st, b_Spower, wt, b_Spower->size[1],
                                    wt->size[0]);
              i1 = Spower_iter->size[0];
              Spower_iter->size[0] = b_Spower->size[0];
              emxEnsureCapacity_real32_T(&e_st, Spower_iter, i1, &bc_emlrtRTEI);
              b_loop_ub = b_Spower->size[0];
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                Spower_iter->data[i1] = 0.0F;
                c_loop_ub = b_Spower->size[1];
                for (end = 0; end < c_loop_ub; end++) {
                  Spower_iter->data[i1] += b_Spower->data[i1 + b_Spower->size[0]
                    * end] * static_cast<real32_T>(wt->data[end]);
                }
              }
            } else {
              /*  uniform weights */
              i1 = wt->size[0];
              wt->size[0] = num_tapers;
              emxEnsureCapacity_real_T(&d_st, wt, i1, &xb_emlrtRTEI);
              for (i1 = 0; i1 < num_tapers; i1++) {
                wt->data[i1] = 1.0 / static_cast<real_T>(num_tapers);
              }

              e_st.site = &r_emlrtRSI;
              f_st.site = &fe_emlrtRSI;
              b_dynamic_size_checks(&f_st, b_Spower, wt, b_Spower->size[1],
                                    wt->size[0]);
              i1 = Spower_iter->size[0];
              Spower_iter->size[0] = b_Spower->size[0];
              emxEnsureCapacity_real32_T(&e_st, Spower_iter, i1, &ac_emlrtRTEI);
              b_loop_ub = b_Spower->size[0];
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                Spower_iter->data[i1] = 0.0F;
                c_loop_ub = b_Spower->size[1];
                for (end = 0; end < c_loop_ub; end++) {
                  Spower_iter->data[i1] += b_Spower->data[i1 + b_Spower->size[0]
                    * end] * static_cast<real32_T>(wt->data[end]);
                }
              }
            }

            /* Append the spectrum to the spectrogram */
            i1 = n + 1;
            if ((i1 < 1) || (i1 > mt_spectrogram->size[1])) {
              emlrtDynamicBoundsCheckR2012b(i1, 1, mt_spectrogram->size[1],
                &h_emlrtBCI, &d_st);
            }

            end = freq_inds->size[1];
            for (c_loop_ub = 0; c_loop_ub < end; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                i1 = c_loop_ub + 1;
                if ((i1 < 1) || (i1 > Spower_iter->size[0])) {
                  emlrtDynamicBoundsCheckR2012b(i1, 1, Spower_iter->size[0],
                    &i_emlrtBCI, &d_st);
                }
              }
            }

            end = freq_inds->size[1] - 1;
            b_loop_ub = 0;
            for (c_loop_ub = 0; c_loop_ub <= end; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                b_loop_ub++;
              }
            }

            i1 = r3->size[0] * r3->size[1];
            r3->size[0] = 1;
            r3->size[1] = b_loop_ub;
            emxEnsureCapacity_int32_T(&d_st, r3, i1, &cb_emlrtRTEI);
            b_loop_ub = 0;
            for (c_loop_ub = 0; c_loop_ub <= end; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                r3->data[b_loop_ub] = c_loop_ub + 1;
                b_loop_ub++;
              }
            }

            b_loop_ub = r3->size[1];
            i1 = r2->size[0];
            r2->size[0] = r3->size[1];
            emxEnsureCapacity_int32_T(&d_st, r2, i1, &fc_emlrtRTEI);
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              r2->data[i1] = r3->data[i1];
            }

            d_window_start[0] = mt_spectrogram->size[0];
            emlrtSubAssignSizeCheckR2012b(&d_window_start[0], 1, &r2->size[0], 1,
              &g_emlrtECI, &d_st);
            end = freq_inds->size[1] - 1;
            b_loop_ub = 0;
            for (c_loop_ub = 0; c_loop_ub <= end; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                b_loop_ub++;
              }
            }

            i1 = r1->size[0] * r1->size[1];
            r1->size[0] = 1;
            r1->size[1] = b_loop_ub;
            emxEnsureCapacity_int32_T(&d_st, r1, i1, &cb_emlrtRTEI);
            b_loop_ub = 0;
            for (c_loop_ub = 0; c_loop_ub <= end; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                r1->data[b_loop_ub] = c_loop_ub + 1;
                b_loop_ub++;
              }
            }

            b_loop_ub = r1->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              mt_spectrogram->data[i1 + mt_spectrogram->size[0] * n] =
                Spower_iter->data[r1->data[i1] - 1];
            }
          }
        }

        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(&d_st);
        }
      } catch (...) {
        emlrtHadParallelError = true;
      }
    }

    if (!emlrtHadParallelError) {
      emlrtHeapReferenceStackLeaveScope(&d_st, 19);
      emxFree_real_T(&wt);
      emxFree_real32_T(&wk);
      emxFree_real32_T(&b_b);
      emxFree_real32_T(&Spower_iter);
      emxFree_real32_T(&a);
      emxFree_real32_T(&b_Spower);
      emxFree_creal32_T(&fft_data);
      emxFree_real32_T(&c_data_segment);
      emxFree_real_T(&c_y);
      emxFree_boolean_T(&b_data_segment);
      emxFree_boolean_T(&r5);
      emxFree_real32_T(&data_segment);
      emxFree_real32_T(&Spower);
      emxFree_real32_T(&b_y);
      emxFree_real32_T(&y);
      emxFree_real_T(&r4);
      emxFree_int32_T(&r3);
      emxFree_int32_T(&r2);
      emxFree_int32_T(&r1);
    }
  }

  emlrtExitParallelRegion(sp, omp_in_parallel());
  emxFree_boolean_T(&freq_inds);

  /* Compute one-sided PSD spectrum  */
  st.site = &s_emlrtRSI;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(&st, x, i, &gb_emlrtRTEI);
  loop_ub = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < loop_ub; i++) {
    x->data[i] = (sfreqs->data[i] == 0.0);
  }

  emxInit_int32_T(&st, &ii, 2, &wc_emlrtRTEI, true);
  b_st.site = &tf_emlrtRSI;
  eml_find(&b_st, x, ii);
  i = window_start->size[0] * window_start->size[1];
  window_start->size[0] = 1;
  window_start->size[1] = ii->size[1];
  emxEnsureCapacity_real_T(&st, window_start, i, &hb_emlrtRTEI);
  loop_ub = ii->size[0] * ii->size[1];
  for (i = 0; i < loop_ub; i++) {
    window_start->data[i] = ii->data[i];
  }

  st.site = &t_emlrtRSI;
  b = Fs / 2.0;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(&st, x, i, &ib_emlrtRTEI);
  loop_ub = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < loop_ub; i++) {
    x->data[i] = (sfreqs->data[i] == b);
  }

  emxInit_real_T(&st, &d_y, 2, &xc_emlrtRTEI, true);
  b_st.site = &tf_emlrtRSI;
  eml_find(&b_st, x, ii);
  emxFree_boolean_T(&x);
  if (sfreqs->size[1] < 1) {
    d_y->size[0] = 1;
    d_y->size[1] = 0;
  } else {
    i = d_y->size[0] * d_y->size[1];
    d_y->size[0] = 1;
    d_y->size[1] = static_cast<int32_T>(static_cast<real_T>(sfreqs->size[1]) -
      1.0) + 1;
    emxEnsureCapacity_real_T(sp, d_y, i, &jb_emlrtRTEI);
    loop_ub = static_cast<int32_T>(static_cast<real_T>(sfreqs->size[1]) - 1.0);
    for (i = 0; i <= loop_ub; i++) {
      d_y->data[i] = static_cast<real_T>(i) + 1.0;
    }
  }

  emxInit_real_T(sp, &b_window_start, 2, &kb_emlrtRTEI, true);
  st.site = &u_emlrtRSI;
  b_st.site = &wf_emlrtRSI;
  i = b_window_start->size[0] * b_window_start->size[1];
  b_window_start->size[0] = 1;
  b_window_start->size[1] = window_start->size[1] + ii->size[1];
  emxEnsureCapacity_real_T(&b_st, b_window_start, i, &kb_emlrtRTEI);
  loop_ub = window_start->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_window_start->data[i] = window_start->data[i];
  }

  loop_ub = ii->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_window_start->data[i + window_start->size[1]] = ii->data[i];
  }

  emxInit_real_T(&b_st, &b_select, 2, &vc_emlrtRTEI, true);
  emxInit_real32_T(&b_st, &c_b, 2, &mb_emlrtRTEI, true);
  emxInit_int32_T(&b_st, &ia, 1, &sb_emlrtRTEI, true);
  c_st.site = &xf_emlrtRSI;
  do_vectors(&c_st, d_y, b_window_start, b_select, ia, ib_size);
  loop_ub = mt_spectrogram->size[1];
  i = c_b->size[0] * c_b->size[1];
  c_b->size[0] = b_select->size[1];
  c_b->size[1] = mt_spectrogram->size[1];
  emxEnsureCapacity_real32_T(sp, c_b, i, &mb_emlrtRTEI);
  emxFree_real_T(&b_window_start);
  emxFree_real_T(&d_y);
  for (i = 0; i < loop_ub; i++) {
    nz = b_select->size[1];
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < nz; input_sizes_idx_0++) {
      if (b_select->data[input_sizes_idx_0] != static_cast<int32_T>
          (muDoubleScalarFloor(b_select->data[input_sizes_idx_0]))) {
        emlrtIntegerCheckR2012b(b_select->data[input_sizes_idx_0], &b_emlrtDCI,
          sp);
      }

      trueCount = static_cast<int32_T>(b_select->data[input_sizes_idx_0]);
      if ((trueCount < 1) || (trueCount > mt_spectrogram->size[0])) {
        emlrtDynamicBoundsCheckR2012b(trueCount, 1, mt_spectrogram->size[0],
          &e_emlrtBCI, sp);
      }

      c_b->data[input_sizes_idx_0 + c_b->size[0] * i] = mt_spectrogram->data
        [(trueCount + mt_spectrogram->size[0] * i) - 1];
    }
  }

  emxFree_real_T(&b_select);
  loop_ub = c_b->size[0] * c_b->size[1];
  for (i = 0; i < loop_ub; i++) {
    c_b->data[i] *= 2.0F;
  }

  emxInit_real32_T(sp, &varargin_1, 2, &ob_emlrtRTEI, true);
  st.site = &v_emlrtRSI;
  loop_ub = mt_spectrogram->size[1];
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = window_start->size[1];
  varargin_1->size[1] = mt_spectrogram->size[1];
  emxEnsureCapacity_real32_T(&st, varargin_1, i, &ob_emlrtRTEI);
  for (i = 0; i < loop_ub; i++) {
    nz = window_start->size[1];
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < nz; input_sizes_idx_0++) {
      trueCount = static_cast<int32_T>(window_start->data[input_sizes_idx_0]);
      if ((trueCount < 1) || (trueCount > mt_spectrogram->size[0])) {
        emlrtDynamicBoundsCheckR2012b(trueCount, 1, mt_spectrogram->size[0],
          &b_emlrtBCI, &st);
      }

      varargin_1->data[input_sizes_idx_0 + varargin_1->size[0] * i] =
        mt_spectrogram->data[(trueCount + mt_spectrogram->size[0] * i) - 1];
    }
  }

  emxInit_real32_T(&st, &varargin_3, 2, &qb_emlrtRTEI, true);
  loop_ub = mt_spectrogram->size[1];
  i = varargin_3->size[0] * varargin_3->size[1];
  varargin_3->size[0] = ii->size[1];
  varargin_3->size[1] = mt_spectrogram->size[1];
  emxEnsureCapacity_real32_T(&st, varargin_3, i, &qb_emlrtRTEI);
  for (i = 0; i < loop_ub; i++) {
    nz = ii->size[1];
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < nz; input_sizes_idx_0++) {
      if ((ii->data[input_sizes_idx_0] < 1) || (ii->data[input_sizes_idx_0] >
           mt_spectrogram->size[0])) {
        emlrtDynamicBoundsCheckR2012b(ii->data[input_sizes_idx_0], 1,
          mt_spectrogram->size[0], &emlrtBCI, &st);
      }

      varargin_3->data[input_sizes_idx_0 + varargin_3->size[0] * i] =
        mt_spectrogram->data[(ii->data[input_sizes_idx_0] + mt_spectrogram->
        size[0] * i) - 1];
    }
  }

  b_st.site = &kg_emlrtRSI;
  i = ia->size[0];
  ia->size[0] = window_start->size[1];
  emxEnsureCapacity_int32_T(&b_st, ia, i, &sb_emlrtRTEI);
  loop_ub = window_start->size[1];
  for (i = 0; i < loop_ub; i++) {
    ia->data[i] = static_cast<int32_T>(window_start->data[i]);
  }

  if ((ia->size[0] != 0) && (mt_spectrogram->size[1] != 0)) {
    b_i = mt_spectrogram->size[1];
  } else if ((c_b->size[0] != 0) && (c_b->size[1] != 0)) {
    b_i = c_b->size[1];
  } else {
    i = ia->size[0];
    ia->size[0] = ii->size[1];
    emxEnsureCapacity_int32_T(&b_st, ia, i, &ub_emlrtRTEI);
    loop_ub = ii->size[1];
    for (i = 0; i < loop_ub; i++) {
      ia->data[i] = ii->data[i];
    }

    if ((ia->size[0] != 0) && (mt_spectrogram->size[1] != 0)) {
      b_i = mt_spectrogram->size[1];
    } else {
      if (mt_spectrogram->size[1] > 0) {
        b_i = mt_spectrogram->size[1];
      } else {
        b_i = 0;
      }

      if (c_b->size[1] > b_i) {
        b_i = c_b->size[1];
      }

      if (mt_spectrogram->size[1] > b_i) {
        b_i = mt_spectrogram->size[1];
      }
    }
  }

  c_st.site = &lg_emlrtRSI;
  if (mt_spectrogram->size[1] != b_i) {
    i = ia->size[0];
    ia->size[0] = window_start->size[1];
    emxEnsureCapacity_int32_T(&c_st, ia, i, &sb_emlrtRTEI);
    loop_ub = window_start->size[1];
    for (i = 0; i < loop_ub; i++) {
      ia->data[i] = static_cast<int32_T>(window_start->data[i]);
    }

    if ((ia->size[0] != 0) && (mt_spectrogram->size[1] != 0)) {
      emlrtErrorWithMessageIdR2018a(&c_st, &d_emlrtRTEI,
        "MATLAB:catenate:matrixDimensionMismatch",
        "MATLAB:catenate:matrixDimensionMismatch", 0);
    }
  }

  if ((c_b->size[1] != b_i) && ((c_b->size[0] != 0) && (c_b->size[1] != 0))) {
    emlrtErrorWithMessageIdR2018a(&c_st, &d_emlrtRTEI,
      "MATLAB:catenate:matrixDimensionMismatch",
      "MATLAB:catenate:matrixDimensionMismatch", 0);
  }

  if (mt_spectrogram->size[1] != b_i) {
    i = ia->size[0];
    ia->size[0] = ii->size[1];
    emxEnsureCapacity_int32_T(&c_st, ia, i, &ub_emlrtRTEI);
    loop_ub = ii->size[1];
    for (i = 0; i < loop_ub; i++) {
      ia->data[i] = ii->data[i];
    }

    if ((ia->size[0] != 0) && (mt_spectrogram->size[1] != 0)) {
      emlrtErrorWithMessageIdR2018a(&c_st, &d_emlrtRTEI,
        "MATLAB:catenate:matrixDimensionMismatch",
        "MATLAB:catenate:matrixDimensionMismatch", 0);
    }
  }

  empty_non_axis_sizes = (b_i == 0);
  guard1 = false;
  if (empty_non_axis_sizes) {
    guard1 = true;
  } else {
    i = ia->size[0];
    ia->size[0] = window_start->size[1];
    emxEnsureCapacity_int32_T(&b_st, ia, i, &sb_emlrtRTEI);
    loop_ub = window_start->size[1];
    for (i = 0; i < loop_ub; i++) {
      ia->data[i] = static_cast<int32_T>(window_start->data[i]);
    }

    if ((ia->size[0] != 0) && (mt_spectrogram->size[1] != 0)) {
      guard1 = true;
    } else {
      input_sizes_idx_0 = 0;
    }
  }

  if (guard1) {
    i = ia->size[0];
    ia->size[0] = window_start->size[1];
    emxEnsureCapacity_int32_T(&b_st, ia, i, &sb_emlrtRTEI);
    loop_ub = window_start->size[1];
    for (i = 0; i < loop_ub; i++) {
      ia->data[i] = static_cast<int32_T>(window_start->data[i]);
    }

    input_sizes_idx_0 = ia->size[0];
  }

  emxFree_real_T(&window_start);
  if (empty_non_axis_sizes || ((c_b->size[0] != 0) && (c_b->size[1] != 0))) {
    num_tapers = c_b->size[0];
  } else {
    num_tapers = 0;
  }

  guard1 = false;
  if (empty_non_axis_sizes) {
    guard1 = true;
  } else {
    i = ia->size[0];
    ia->size[0] = ii->size[1];
    emxEnsureCapacity_int32_T(&b_st, ia, i, &ub_emlrtRTEI);
    loop_ub = ii->size[1];
    for (i = 0; i < loop_ub; i++) {
      ia->data[i] = ii->data[i];
    }

    if ((ia->size[0] != 0) && (mt_spectrogram->size[1] != 0)) {
      guard1 = true;
    } else {
      nz = 0;
    }
  }

  if (guard1) {
    i = ia->size[0];
    ia->size[0] = ii->size[1];
    emxEnsureCapacity_int32_T(&b_st, ia, i, &ub_emlrtRTEI);
    loop_ub = ii->size[1];
    for (i = 0; i < loop_ub; i++) {
      ia->data[i] = ii->data[i];
    }

    nz = ia->size[0];
  }

  emxFree_int32_T(&ia);
  emxFree_int32_T(&ii);
  trueCount = input_sizes_idx_0;
  i = mt_spectrogram->size[0] * mt_spectrogram->size[1];
  mt_spectrogram->size[0] = b_i;
  mt_spectrogram->size[1] = (trueCount + num_tapers) + nz;
  emxEnsureCapacity_real32_T(sp, mt_spectrogram, i, &dc_emlrtRTEI);
  for (i = 0; i < trueCount; i++) {
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < b_i; input_sizes_idx_0++) {
      mt_spectrogram->data[input_sizes_idx_0 + mt_spectrogram->size[0] * i] =
        varargin_1->data[i + trueCount * input_sizes_idx_0] /
        static_cast<real32_T>(Fs);
    }
  }

  emxFree_real32_T(&varargin_1);
  for (i = 0; i < num_tapers; i++) {
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < b_i; input_sizes_idx_0++) {
      mt_spectrogram->data[input_sizes_idx_0 + mt_spectrogram->size[0] * (i +
        trueCount)] = c_b->data[i + num_tapers * input_sizes_idx_0] /
        static_cast<real32_T>(Fs);
    }
  }

  emxFree_real32_T(&c_b);
  for (i = 0; i < nz; i++) {
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < b_i; input_sizes_idx_0++) {
      mt_spectrogram->data[input_sizes_idx_0 + mt_spectrogram->size[0] * ((i +
        trueCount) + num_tapers)] = varargin_3->data[i + nz * input_sizes_idx_0]
        / static_cast<real32_T>(Fs);
    }
  }

  emxFree_real32_T(&varargin_3);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (multitaper_spectrogram_coder.cpp) */
