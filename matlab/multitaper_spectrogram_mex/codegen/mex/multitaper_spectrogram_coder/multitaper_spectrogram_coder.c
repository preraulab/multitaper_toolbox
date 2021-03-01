/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder.c
 *
 * Code generation for function 'multitaper_spectrogram_coder'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder.h"
#include "all.h"
#include "any.h"
#include "colon.h"
#include "combineVectorElements.h"
#include "detrend.h"
#include "eml_mtimes_helper.h"
#include "eml_setop.h"
#include "fft.h"
#include "find.h"
#include "indexShapeCheck.h"
#include "mean.h"
#include "mtimes.h"
#include "multitaper_spectrogram_coder_data.h"
#include "multitaper_spectrogram_coder_emxutil.h"
#include "multitaper_spectrogram_coder_types.h"
#include "nextpow2.h"
#include "power.h"
#include "repmat.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 55,    /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 61,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 65,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 77,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 98,  /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 101, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 105, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 112, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 114, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 118, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 121, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 124, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 129, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 130, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 134, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 136, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 137, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 142, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 150, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 151, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 152, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 153, /* lineNo */
  "multitaper_spectrogram_coder",      /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 103, /* lineNo */
  "colon",                             /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 15, /* lineNo */
  "sum",                               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

static emlrtRSInfo gf_emlrtRSI = { 39, /* lineNo */
  "find",                              /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pathName */
};

static emlrtRSInfo jf_emlrtRSI = { 19, /* lineNo */
  "setdiff",                           /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/setdiff.m"/* pathName */
};

static emlrtRSInfo kf_emlrtRSI = { 70, /* lineNo */
  "eml_setop",                         /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/private/eml_setop.m"/* pathName */
};

static emlrtRSInfo wf_emlrtRSI = { 27, /* lineNo */
  "cat",                               /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/cat.m"/* pathName */
};

static emlrtRSInfo xf_emlrtRSI = { 102,/* lineNo */
  "cat_impl",                          /* fcnName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/cat.m"/* pathName */
};

static emlrtRTEInfo b_emlrtRTEI = { 283,/* lineNo */
  27,                                  /* colNo */
  "check_non_axis_size",               /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/cat.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 20,/* lineNo */
  15,                                  /* colNo */
  "rdivide_helper",                    /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/eml/+coder/+internal/rdivide_helper.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  146,                                 /* lineNo */
  39,                                  /* colNo */
  "mt_spectrum",                       /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo emlrtECI = { -1,    /* nDims */
  146,                                 /* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  146,                                 /* lineNo */
  22,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  106,                                 /* lineNo */
  26,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo b_emlrtECI = { 2,   /* nDims */
  137,                                 /* lineNo */
  29,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtECInfo c_emlrtECI = { 2,   /* nDims */
  136,                                 /* lineNo */
  16,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtDCInfo emlrtDCI = { 136,   /* lineNo */
  30,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  1                                    /* checkKind */
};

static emlrtECInfo d_emlrtECI = { 2,   /* nDims */
  134,                                 /* lineNo */
  50,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtDCInfo b_emlrtDCI = { 134, /* lineNo */
  79,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  1                                    /* checkKind */
};

static emlrtECInfo e_emlrtECI = { 2,   /* nDims */
  124,                                 /* lineNo */
  14,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtECInfo f_emlrtECI = { 2,   /* nDims */
  118,                                 /* lineNo */
  20,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtECInfo g_emlrtECI = { 2,   /* nDims */
  68,                                  /* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtDCInfo c_emlrtDCI = { 77,  /* lineNo */
  24,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  69,                                  /* lineNo */
  10,                                  /* colNo */
  "sfreqs",                            /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  98,                                  /* lineNo */
  25,                                  /* colNo */
  "window_start",                      /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  98,                                  /* lineNo */
  25,                                  /* colNo */
  "data",                              /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = { 98,  /* lineNo */
  25,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  153,                                 /* lineNo */
  65,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = { 153, /* lineNo */
  65,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  153,                                 /* lineNo */
  34,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  153,                                 /* lineNo */
  91,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  130,                                 /* lineNo */
  37,                                  /* colNo */
  "Spower",                            /* aName */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo x_emlrtRTEI = { 55,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo y_emlrtRTEI = { 65,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = { 68,/* lineNo */
  14,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = { 68,/* lineNo */
  47,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = { 68,/* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo db_emlrtRTEI = { 1,/* lineNo */
  45,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = { 73,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = { 77,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = { 88,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = { 92,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = { 90,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo jb_emlrtRTEI = { 150,/* lineNo */
  18,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo kb_emlrtRTEI = { 150,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = { 151,/* lineNo */
  23,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo mb_emlrtRTEI = { 98,/* lineNo */
  44,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo nb_emlrtRTEI = { 28,/* lineNo */
  9,                                   /* colNo */
  "colon",                             /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/colon.m"/* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = { 98,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = { 152,/* lineNo */
  36,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = { 101,/* lineNo */
  12,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = { 15,/* lineNo */
  13,                                  /* colNo */
  "isnan",                             /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/isnan.m"/* pName */
};

static emlrtRTEInfo sb_emlrtRTEI = { 153,/* lineNo */
  50,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo tb_emlrtRTEI = { 112,/* lineNo */
  32,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo ub_emlrtRTEI = { 153,/* lineNo */
  19,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = { 153,/* lineNo */
  76,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = { 124,/* lineNo */
  14,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = { 124,/* lineNo */
  34,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo yb_emlrtRTEI = { 142,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo ac_emlrtRTEI = { 130,/* lineNo */
  28,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo bc_emlrtRTEI = { 131,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo cc_emlrtRTEI = { 153,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo dc_emlrtRTEI = { 134,/* lineNo */
  50,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo ec_emlrtRTEI = { 134,/* lineNo */
  74,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo fc_emlrtRTEI = { 76,/* lineNo */
  9,                                   /* colNo */
  "eml_mtimes_helper",                 /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pName */
};

static emlrtRTEInfo gc_emlrtRTEI = { 136,/* lineNo */
  25,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo hc_emlrtRTEI = { 137,/* lineNo */
  29,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo ic_emlrtRTEI = { 137,/* lineNo */
  25,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo jc_emlrtRTEI = { 136,/* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo kc_emlrtRTEI = { 134,/* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo lc_emlrtRTEI = { 130,/* lineNo */
  9,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo mc_emlrtRTEI = { 124,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo nc_emlrtRTEI = { 121,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo oc_emlrtRTEI = { 137,/* lineNo */
  47,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo pc_emlrtRTEI = { 146,/* lineNo */
  39,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo qc_emlrtRTEI = { 68,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo rc_emlrtRTEI = { 152,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

static emlrtRTEInfo sc_emlrtRTEI = { 31,/* lineNo */
  6,                                   /* colNo */
  "find",                              /* fName */
  "/apps/software/MATLAB/R2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pName */
};

static emlrtRTEInfo tc_emlrtRTEI = { 152,/* lineNo */
  18,                                  /* colNo */
  "multitaper_spectrogram_coder",      /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder.m"/* pName */
};

/* Function Definitions */
emlrtCTX emlrtGetRootTLSGlobal(void)
{
  return emlrtRootTLSGlobal;
}

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, const emlrtConstCTX aTLS,
  void *aData)
{
  omp_set_lock(&emlrtLockGlobal);
  emlrtCallLockeeFunction(aLockee, aTLS, aData);
  omp_unset_lock(&emlrtLockGlobal);
}

void multitaper_spectrogram_coder(const emlrtStack *sp, const emxArray_real32_T *
  data, real_T Fs, const real_T frequency_range[2], const emxArray_real_T
  *DPSS_tapers, const emxArray_real_T *DPSS_eigen, real_T winstep_samples,
  real_T min_NFFT, real_T detrend_opt, real_T weighting, emxArray_real32_T
  *mt_spectrogram, emxArray_real_T *stimes, emxArray_real_T *sfreqs)
{
  jmp_buf b_emlrtJBEnviron;
  jmp_buf * volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_boolean_T *Spower_iter;
  emxArray_boolean_T *b_sfreqs;
  emxArray_boolean_T *freq_inds;
  emxArray_boolean_T *r;
  emxArray_creal32_T *fft_data;
  emxArray_int32_T *ia;
  emxArray_int32_T *ii;
  emxArray_int32_T *r1;
  emxArray_int32_T *r2;
  emxArray_real32_T *Spower;
  emxArray_real32_T *a;
  emxArray_real32_T *b_Spower;
  emxArray_real32_T *b_Spower_iter;
  emxArray_real32_T *b_b;
  emxArray_real32_T *b_y;
  emxArray_real32_T *c_b;
  emxArray_real32_T *r3;
  emxArray_real32_T *varargin_1;
  emxArray_real32_T *varargin_3;
  emxArray_real32_T *wk;
  emxArray_real32_T *y;
  emxArray_real_T *b_select;
  emxArray_real_T *b_window_start;
  emxArray_real_T *c_y;
  emxArray_real_T *r4;
  emxArray_real_T *r5;
  emxArray_real_T *window_start;
  emxArray_real_T *wt;
  real_T b;
  real_T c_window_start;
  real_T d;
  real_T nfft;
  int32_T b_wk[2];
  int32_T c_Spower[2];
  int32_T ib_size[1];
  int32_T b_i;
  int32_T b_ii;
  int32_T b_input_sizes_idx_0;
  int32_T b_loop_ub;
  int32_T c_loop_ub;
  int32_T d_loop_ub;
  int32_T e_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T input_sizes_idx_0;
  int32_T loop_ub;
  int32_T n;
  int32_T num_tapers;
  int32_T nz;
  int32_T winsize_samples;
  real32_T Tpower;
  boolean_T emlrtHadParallelError = false;
  boolean_T empty_non_axis_sizes;
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
  /*    data: <number of samples> x 1 vector - time series data -- required */
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
  /*    mt_spectrogram: FxT matrix of one-sided power spectral density (PSD) */
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
  winsize_samples = DPSS_tapers->size[0];
  num_tapers = DPSS_tapers->size[1];

  /* Total data length */
  /* Window start indices */
  st.site = &emlrtRSI;
  b = (real_T)(data->size[0] - DPSS_tapers->size[0]) + 1.0;
  emxInit_real_T(&st, &window_start, 2, &x_emlrtRTEI, true);
  if (muDoubleScalarIsNaN(winstep_samples)) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &x_emlrtRTEI);
    window_start->data[0] = rtNaN;
  } else if ((winstep_samples == 0.0) || ((1.0 < b) && (winstep_samples < 0.0)) ||
             ((b < 1.0) && (winstep_samples > 0.0))) {
    window_start->size[0] = 1;
    window_start->size[1] = 0;
  } else if (muDoubleScalarIsInf(winstep_samples)) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &x_emlrtRTEI);
    window_start->data[0] = 1.0;
  } else if (muDoubleScalarFloor(winstep_samples) == winstep_samples) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    input_sizes_idx_0 = (int32_T)muDoubleScalarFloor((b - 1.0) / winstep_samples);
    window_start->size[1] = input_sizes_idx_0 + 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &x_emlrtRTEI);
    for (i = 0; i <= input_sizes_idx_0; i++) {
      window_start->data[i] = winstep_samples * (real_T)i + 1.0;
    }
  } else {
    b_st.site = &w_emlrtRSI;
    eml_float_colon(&b_st, winstep_samples, b, window_start);
  }

  /* Number of windows */
  /* Number of points in the FFT */
  st.site = &b_emlrtRSI;
  b_st.site = &b_emlrtRSI;
  b = nextpow2(DPSS_tapers->size[0]);
  b_st.site = &gb_emlrtRSI;
  st.site = &b_emlrtRSI;
  b_st.site = &b_emlrtRSI;
  nfft = nextpow2(min_NFFT);
  b_st.site = &gb_emlrtRSI;
  nfft = muDoubleScalarMax(muDoubleScalarMax(muDoubleScalarPower(2.0, b),
    DPSS_tapers->size[0]), muDoubleScalarPower(2.0, nfft));

  /* Create the frequency vector */
  b = Fs / nfft;
  st.site = &c_emlrtRSI;
  if (muDoubleScalarIsNaN(b) || muDoubleScalarIsNaN(Fs)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &y_emlrtRTEI);
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
    emxEnsureCapacity_real_T(&st, sfreqs, i, &y_emlrtRTEI);
    sfreqs->data[0] = rtNaN;
  } else if (muDoubleScalarIsInf(b)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &y_emlrtRTEI);
    sfreqs->data[0] = 0.0;
  } else if (muDoubleScalarFloor(b) == b) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    input_sizes_idx_0 = (int32_T)muDoubleScalarFloor(Fs / b);
    sfreqs->size[1] = input_sizes_idx_0 + 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &y_emlrtRTEI);
    for (i = 0; i <= input_sizes_idx_0; i++) {
      sfreqs->data[i] = b * (real_T)i;
    }
  } else {
    b_st.site = &w_emlrtRSI;
    b_eml_float_colon(&b_st, b, Fs, sfreqs);
  }

  emxInit_boolean_T(&st, &b_sfreqs, 2, &jb_emlrtRTEI, true);

  /*  all possible frequencies */
  /* Get just the frequencies for the given frequency range */
  i = b_sfreqs->size[0] * b_sfreqs->size[1];
  b_sfreqs->size[0] = 1;
  b_sfreqs->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, b_sfreqs, i, &ab_emlrtRTEI);
  input_sizes_idx_0 = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    b_sfreqs->data[i] = (sfreqs->data[i] >= frequency_range[0]);
  }

  emxInit_boolean_T(sp, &r, 2, &db_emlrtRTEI, true);
  i = r->size[0] * r->size[1];
  r->size[0] = 1;
  r->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, r, i, &bb_emlrtRTEI);
  input_sizes_idx_0 = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    r->data[i] = (sfreqs->data[i] <= frequency_range[1]);
  }

  emxInit_boolean_T(sp, &freq_inds, 2, &qc_emlrtRTEI, true);
  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])b_sfreqs->size, *(int32_T (*)[2])
    r->size, &g_emlrtECI, sp);
  i = freq_inds->size[0] * freq_inds->size[1];
  freq_inds->size[0] = 1;
  freq_inds->size[1] = b_sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, freq_inds, i, &cb_emlrtRTEI);
  input_sizes_idx_0 = b_sfreqs->size[0] * b_sfreqs->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    freq_inds->data[i] = (b_sfreqs->data[i] && r->data[i]);
  }

  nz = b_sfreqs->size[1] - 1;
  input_sizes_idx_0 = 0;
  for (b_i = 0; b_i <= nz; b_i++) {
    if (b_sfreqs->data[b_i] && r->data[b_i]) {
      input_sizes_idx_0++;
    }
  }

  b_input_sizes_idx_0 = 0;
  for (b_i = 0; b_i <= nz; b_i++) {
    if (b_sfreqs->data[b_i] && r->data[b_i]) {
      if ((b_i + 1 < 1) || (b_i + 1 > sfreqs->size[1])) {
        emlrtDynamicBoundsCheckR2012b(b_i + 1, 1, sfreqs->size[1], &d_emlrtBCI,
          sp);
      }

      sfreqs->data[b_input_sizes_idx_0] = sfreqs->data[b_i];
      b_input_sizes_idx_0++;
    }
  }

  emxFree_boolean_T(&r);
  i = sfreqs->size[0] * sfreqs->size[1];
  sfreqs->size[0] = 1;
  sfreqs->size[1] = input_sizes_idx_0;
  emxEnsureCapacity_real_T(sp, sfreqs, i, &db_emlrtRTEI);

  /* Compute the times of the middle of each spectrum */
  nz = (int32_T)muDoubleScalarRound((real_T)DPSS_tapers->size[0] / 2.0);
  i = stimes->size[0] * stimes->size[1];
  stimes->size[0] = 1;
  stimes->size[1] = window_start->size[1];
  emxEnsureCapacity_real_T(sp, stimes, i, &eb_emlrtRTEI);
  input_sizes_idx_0 = window_start->size[0] * window_start->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    stimes->data[i] = ((window_start->data[i] + (real_T)nz) - 1.0) / Fs;
  }

  /*  stimes start from 0 */
  /* Preallocate spectrogram and slice data for efficient parallel computing */
  st.site = &d_emlrtRSI;
  b_st.site = &ib_emlrtRSI;
  c_st.site = &jb_emlrtRSI;
  nz = combineVectorElements(&c_st, freq_inds);
  if (nz < 0) {
    emlrtNonNegativeCheckR2012b(nz, &c_emlrtDCI, sp);
  }

  i = mt_spectrogram->size[0] * mt_spectrogram->size[1];
  mt_spectrogram->size[0] = nz;
  mt_spectrogram->size[1] = window_start->size[1];
  emxEnsureCapacity_real32_T(sp, mt_spectrogram, i, &fb_emlrtRTEI);
  input_sizes_idx_0 = nz * window_start->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    mt_spectrogram->data[i] = 0.0F;
  }

  /*  COMPUTE THE MULTITAPER SPECTROGRAM */
  /*  */
  /*      STEP 1: Compute DPSS tapers based on desired spectral properties */
  /*      STEP 2: Multiply the data segment by the DPSS Tapers */
  /*      STEP 3: Compute the spectrum for each tapered segment */
  /*      STEP 4: Take the mean of the tapered spectra */
  /*  pre-compute weights */
  emxInit_real_T(sp, &wt, 1, &gb_emlrtRTEI, true);
  if (weighting == 1.0) {
    i = wt->size[0];
    wt->size[0] = DPSS_eigen->size[0];
    emxEnsureCapacity_real_T(sp, wt, i, &gb_emlrtRTEI);
    input_sizes_idx_0 = DPSS_eigen->size[0];
    for (i = 0; i < input_sizes_idx_0; i++) {
      wt->data[i] = DPSS_eigen->data[i] / (real_T)DPSS_tapers->size[1];
    }
  } else if (weighting == 0.0) {
    i = wt->size[0];
    wt->size[0] = DPSS_tapers->size[1];
    emxEnsureCapacity_real_T(sp, wt, i, &ib_emlrtRTEI);
    input_sizes_idx_0 = DPSS_tapers->size[1];
    for (i = 0; i < input_sizes_idx_0; i++) {
      wt->data[i] = 1.0 / (real_T)DPSS_tapers->size[1];
    }
  } else {
    i = wt->size[0];
    wt->size[0] = 1;
    emxEnsureCapacity_real_T(sp, wt, i, &hb_emlrtRTEI);
    wt->data[0] = 0.0;
  }

  nz = window_start->size[1] - 1;
  emlrtEnterParallelRegion(sp, omp_in_parallel());
  emlrtPushJmpBuf(sp, &emlrtJBStack);

#pragma omp parallel \
 num_threads(emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs())) \
 private(r1,r2,y,r3,r4,b_y,Spower,Spower_iter,r5,fft_data,b_Spower,a,b_Spower_iter,b_b,wk,Tpower,b_emlrtJBEnviron,f_st,i1,n,c_window_start,loop_ub,d,b_loop_ub,c_Spower,c_loop_ub,b_ii,i2,d_loop_ub,i3,b_wk,e_loop_ub) \
 firstprivate(d_st,e_st,emlrtHadParallelError)

  {
    if (setjmp(b_emlrtJBEnviron) == 0) {
      d_st.prev = sp;
      d_st.tls = emlrtAllocTLS(sp, omp_get_thread_num());
      d_st.site = NULL;
      emlrtSetJmpBuf(&d_st, &b_emlrtJBEnviron);
      e_st.prev = &d_st;
      e_st.tls = d_st.tls;
      f_st.prev = &e_st;
      f_st.tls = e_st.tls;
      emxInit_int32_T(&d_st, &r1, 2, &pc_emlrtRTEI, true);
      emxInit_int32_T(&d_st, &r2, 2, &pc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &y, 1, &oc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &r3, 2, &ic_emlrtRTEI, true);
      emxInit_real_T(&d_st, &r4, 2, &db_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &b_y, 2, &dc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &Spower, 2, &ac_emlrtRTEI, true);
      emxInit_boolean_T(&d_st, &Spower_iter, 1, &qb_emlrtRTEI, true);
      emxInit_real_T(&d_st, &r5, 2, &db_emlrtRTEI, true);
      emxInit_creal32_T(&d_st, &fft_data, 2, &nc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &b_Spower, 2, &mc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &a, 1, &bc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &b_Spower_iter, 1, &lc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &b_b, 2, &kc_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &wk, 2, &jc_emlrtRTEI, true);
    } else {
      emlrtHadParallelError = true;
    }

#pragma omp for nowait

    for (n = 0; n <= nz; n++) {
      if (emlrtHadParallelError)
        continue;
      if (setjmp(b_emlrtJBEnviron) == 0) {
        /* Loop in parallel over all of the windows */
        /* Grab the data for the given window */
        if (winsize_samples - 1 < 0) {
          r5->size[0] = 1;
          r5->size[1] = 0;
        } else {
          i1 = r5->size[0] * r5->size[1];
          r5->size[0] = 1;
          r5->size[1] = (int32_T)((real_T)winsize_samples - 1.0) + 1;
          emxEnsureCapacity_real_T(&d_st, r5, i1, &mb_emlrtRTEI);
          loop_ub = (int32_T)((real_T)winsize_samples - 1.0);
          for (i1 = 0; i1 <= loop_ub; i1++) {
            r5->data[i1] = i1;
          }
        }

        e_st.site = &e_emlrtRSI;
        indexShapeCheck(&e_st, data->size[0], *(int32_T (*)[2])r5->size);
        if ((n + 1 < 1) || (n + 1 > window_start->size[1])) {
          emlrtDynamicBoundsCheckR2012b(n + 1, 1, window_start->size[1],
            &e_emlrtBCI, &d_st);
        }

        c_window_start = window_start->data[n];
        i1 = b_Spower_iter->size[0];
        b_Spower_iter->size[0] = r5->size[1];
        emxEnsureCapacity_real32_T(&d_st, b_Spower_iter, i1, &ob_emlrtRTEI);
        loop_ub = r5->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          d = c_window_start + r5->data[i1];
          if (d != (int32_T)muDoubleScalarFloor(d)) {
            emlrtIntegerCheckR2012b(d, &d_emlrtDCI, &d_st);
          }

          if (((int32_T)d < 1) || ((int32_T)d > data->size[0])) {
            emlrtDynamicBoundsCheckR2012b((int32_T)d, 1, data->size[0],
              &f_emlrtBCI, &d_st);
          }

          b_Spower_iter->data[i1] = data->data[(int32_T)d - 1];
        }

        /* Skip empty segments */
        i1 = Spower_iter->size[0];
        Spower_iter->size[0] = b_Spower_iter->size[0];
        emxEnsureCapacity_boolean_T(&d_st, Spower_iter, i1, &qb_emlrtRTEI);
        loop_ub = b_Spower_iter->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          Spower_iter->data[i1] = (b_Spower_iter->data[i1] == 0.0F);
        }

        e_st.site = &f_emlrtRSI;
        if (!all(&e_st, Spower_iter)) {
          i1 = Spower_iter->size[0];
          Spower_iter->size[0] = b_Spower_iter->size[0];
          emxEnsureCapacity_boolean_T(&d_st, Spower_iter, i1, &rb_emlrtRTEI);
          loop_ub = b_Spower_iter->size[0];
          for (i1 = 0; i1 < loop_ub; i1++) {
            Spower_iter->data[i1] = muSingleScalarIsNaN(b_Spower_iter->data[i1]);
          }

          e_st.site = &g_emlrtRSI;
          if (any(&e_st, Spower_iter)) {
            if ((n + 1 < 1) || (n + 1 > mt_spectrogram->size[1])) {
              emlrtDynamicBoundsCheckR2012b(n + 1, 1, mt_spectrogram->size[1],
                &c_emlrtBCI, &d_st);
            }

            loop_ub = mt_spectrogram->size[0];
            for (i1 = 0; i1 < loop_ub; i1++) {
              mt_spectrogram->data[i1 + mt_spectrogram->size[0] * n] = rtNaNF;
            }
          } else {
            /* Option to detrend_opt data to remove low frequency DC component */
            if (detrend_opt == 1.0) {
              i1 = a->size[0];
              a->size[0] = b_Spower_iter->size[0];
              emxEnsureCapacity_real32_T(&d_st, a, i1, &tb_emlrtRTEI);
              loop_ub = b_Spower_iter->size[0] - 1;
              for (i1 = 0; i1 <= loop_ub; i1++) {
                a->data[i1] = b_Spower_iter->data[i1];
              }

              e_st.site = &h_emlrtRSI;
              detrend(&e_st, a, b_Spower_iter);
            } else {
              if (detrend_opt == 2.0) {
                e_st.site = &i_emlrtRSI;
                b_detrend(&e_st, b_Spower_iter);
              }
            }

            /* Multiply the data by the tapers (STEP 2) */
            e_st.site = &j_emlrtRSI;
            repmat(&e_st, b_Spower_iter, num_tapers, wk);
            emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])wk->size, *(int32_T (*)[2])
              DPSS_tapers->size, &f_emlrtECI, &d_st);
            loop_ub = wk->size[0] * wk->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              wk->data[i1] *= (real32_T)DPSS_tapers->data[i1];
            }

            /* Compute the FFT (STEP 3) */
            e_st.site = &k_emlrtRSI;
            fft(&e_st, wk, nfft, fft_data);

            /* Compute the weighted mean spectral power across tapers (STEP 4) */
            i1 = b_b->size[0] * b_b->size[1];
            b_b->size[0] = fft_data->size[0];
            b_b->size[1] = fft_data->size[1];
            emxEnsureCapacity_real32_T(&d_st, b_b, i1, &wb_emlrtRTEI);
            loop_ub = fft_data->size[0] * fft_data->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              b_b->data[i1] = fft_data->data[i1].im;
            }

            e_st.site = &l_emlrtRSI;
            power(&e_st, b_b, b_Spower);
            i1 = b_b->size[0] * b_b->size[1];
            b_b->size[0] = fft_data->size[0];
            b_b->size[1] = fft_data->size[1];
            emxEnsureCapacity_real32_T(&d_st, b_b, i1, &xb_emlrtRTEI);
            loop_ub = fft_data->size[0] * fft_data->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              b_b->data[i1] = fft_data->data[i1].re;
            }

            e_st.site = &l_emlrtRSI;
            power(&e_st, b_b, wk);
            emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])b_Spower->size, *(int32_T
              (*)[2])wk->size, &e_emlrtECI, &d_st);
            loop_ub = b_Spower->size[0] * b_Spower->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              b_Spower->data[i1] += wk->data[i1];
            }

            if (weighting == 2.0) {
              /*  adaptive weights - for colored noise spectrum (Percival & Walden */
              /*  p368-p370) */
              e_st.site = &m_emlrtRSI;
              f_st.site = &td_emlrtRSI;
              dynamic_size_checks(&f_st, b_Spower_iter, b_Spower_iter,
                                  b_Spower_iter->size[0], b_Spower_iter->size[0]);
              Tpower = mtimes(b_Spower_iter, b_Spower_iter) / (real32_T)
                b_Spower_iter->size[0];
              loop_ub = b_Spower->size[0];
              i1 = Spower->size[0] * Spower->size[1];
              Spower->size[0] = b_Spower->size[0];
              Spower->size[1] = 2;
              emxEnsureCapacity_real32_T(&d_st, Spower, i1, &ac_emlrtRTEI);
              for (i1 = 0; i1 < loop_ub; i1++) {
                if (1 > b_Spower->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(1, 1, b_Spower->size[1],
                    &j_emlrtBCI, &d_st);
                }

                Spower->data[i1] = b_Spower->data[i1];
              }

              for (i1 = 0; i1 < loop_ub; i1++) {
                if (2 > b_Spower->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(2, 1, b_Spower->size[1],
                    &j_emlrtBCI, &d_st);
                }

                Spower->data[i1 + Spower->size[0]] = b_Spower->data[i1 +
                  b_Spower->size[0]];
              }

              e_st.site = &n_emlrtRSI;
              mean(&e_st, Spower, b_Spower_iter);
              loop_ub = DPSS_eigen->size[0];
              i1 = a->size[0];
              a->size[0] = DPSS_eigen->size[0];
              emxEnsureCapacity_real32_T(&d_st, a, i1, &bc_emlrtRTEI);
              for (i1 = 0; i1 < loop_ub; i1++) {
                a->data[i1] = (real32_T)(1.0 - DPSS_eigen->data[i1]) * Tpower;
              }

              loop_ub = DPSS_eigen->size[0];
              b_loop_ub = a->size[0];
              c_Spower[0] = b_Spower->size[1];
              c_Spower[1] = b_Spower->size[0];
              i1 = (int32_T)muDoubleScalarFloor(nfft);
              c_loop_ub = (int32_T)nfft;
              for (b_ii = 0; b_ii < 3; b_ii++) {
                /*  run 3 iterations */
                /*  calculate the MSE weights */
                i2 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = b_Spower_iter->size[0];
                b_y->size[1] = DPSS_eigen->size[0];
                emxEnsureCapacity_real32_T(&d_st, b_y, i2, &dc_emlrtRTEI);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  d_loop_ub = b_Spower_iter->size[0];
                  for (i3 = 0; i3 < d_loop_ub; i3++) {
                    b_y->data[i3 + b_y->size[0] * i2] = b_Spower_iter->data[i3] *
                      (real32_T)DPSS_eigen->data[i2];
                  }
                }

                if (nfft != i1) {
                  emlrtIntegerCheckR2012b(nfft, &b_emlrtDCI, &d_st);
                }

                i2 = wk->size[0] * wk->size[1];
                wk->size[0] = c_loop_ub;
                wk->size[1] = a->size[0];
                emxEnsureCapacity_real32_T(&d_st, wk, i2, &ec_emlrtRTEI);
                for (i2 = 0; i2 < b_loop_ub; i2++) {
                  for (i3 = 0; i3 < c_loop_ub; i3++) {
                    wk->data[i3 + wk->size[0] * i2] = a->data[i2];
                  }
                }

                emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])b_y->size, *(int32_T (*)
                  [2])wk->size, &d_emlrtECI, &d_st);
                i2 = b_b->size[0] * b_b->size[1];
                b_b->size[0] = b_Spower_iter->size[0];
                b_b->size[1] = num_tapers;
                emxEnsureCapacity_real32_T(&d_st, b_b, i2, &fc_emlrtRTEI);
                for (i2 = 0; i2 < num_tapers; i2++) {
                  d_loop_ub = b_Spower_iter->size[0];
                  for (i3 = 0; i3 < d_loop_ub; i3++) {
                    b_b->data[i3 + b_b->size[0] * i2] = b_Spower_iter->data[i3];
                  }
                }

                e_st.site = &o_emlrtRSI;
                d_loop_ub = b_y->size[0] * b_y->size[1];
                for (i2 = 0; i2 < d_loop_ub; i2++) {
                  b_y->data[i2] += wk->data[i2];
                }

                if ((b_b->size[0] != b_y->size[0]) || (b_b->size[1] != b_y->
                     size[1])) {
                  emlrtErrorWithMessageIdR2018a(&e_st, &d_emlrtRTEI,
                    "MATLAB:dimagree", "MATLAB:dimagree", 0);
                }

                d_loop_ub = b_b->size[0] * b_b->size[1];
                for (i2 = 0; i2 < d_loop_ub; i2++) {
                  b_b->data[i2] /= b_y->data[i2];
                }

                /*  calculate new spectral estimate */
                e_st.site = &p_emlrtRSI;
                power(&e_st, b_b, wk);
                if (c_loop_ub != i1) {
                  emlrtIntegerCheckR2012b(nfft, &emlrtDCI, &d_st);
                }

                i2 = r4->size[0] * r4->size[1];
                r4->size[0] = c_loop_ub;
                r4->size[1] = DPSS_eigen->size[0];
                emxEnsureCapacity_real_T(&d_st, r4, i2, &gc_emlrtRTEI);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  for (i3 = 0; i3 < c_loop_ub; i3++) {
                    r4->data[i3 + r4->size[0] * i2] = DPSS_eigen->data[i2];
                  }
                }

                emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])wk->size, *(int32_T (*)
                  [2])r4->size, &c_emlrtECI, &d_st);
                d_loop_ub = wk->size[0] * wk->size[1];
                for (i2 = 0; i2 < d_loop_ub; i2++) {
                  wk->data[i2] *= (real32_T)r4->data[i2];
                }

                b_wk[0] = wk->size[1];
                b_wk[1] = wk->size[0];
                if ((wk->size[1] != c_Spower[0]) || (wk->size[0] != c_Spower[1]))
                {
                  emlrtSizeEqCheckNDR2012b(&b_wk[0], &c_Spower[0], &b_emlrtECI,
                    &d_st);
                }

                e_st.site = &q_emlrtRSI;
                i2 = b_b->size[0] * b_b->size[1];
                b_b->size[0] = wk->size[1];
                b_b->size[1] = wk->size[0];
                emxEnsureCapacity_real32_T(&e_st, b_b, i2, &hc_emlrtRTEI);
                d_loop_ub = wk->size[0];
                for (i2 = 0; i2 < d_loop_ub; i2++) {
                  e_loop_ub = wk->size[1];
                  for (i3 = 0; i3 < e_loop_ub; i3++) {
                    b_b->data[i3 + b_b->size[0] * i2] = wk->data[i2 + wk->size[0]
                      * i3] * b_Spower->data[i2 + b_Spower->size[0] * i3];
                  }
                }

                f_st.site = &q_emlrtRSI;
                sum(&f_st, b_b, r3);
                i2 = b_Spower_iter->size[0];
                b_Spower_iter->size[0] = r3->size[1];
                emxEnsureCapacity_real32_T(&e_st, b_Spower_iter, i2,
                  &ic_emlrtRTEI);
                d_loop_ub = r3->size[1];
                for (i2 = 0; i2 < d_loop_ub; i2++) {
                  b_Spower_iter->data[i2] = r3->data[i2];
                }

                f_st.site = &q_emlrtRSI;
                b_sum(&f_st, wk, y);
                d_loop_ub = b_Spower_iter->size[0];
                if (b_Spower_iter->size[0] != y->size[0]) {
                  emlrtErrorWithMessageIdR2018a(&e_st, &d_emlrtRTEI,
                    "MATLAB:dimagree", "MATLAB:dimagree", 0);
                }

                for (i2 = 0; i2 < d_loop_ub; i2++) {
                  b_Spower_iter->data[i2] /= y->data[i2];
                }

                if (*emlrtBreakCheckR2012bFlagVar != 0) {
                  emlrtBreakCheckR2012b(&d_st);
                }
              }
            } else {
              /*  eigenvalue or uniform weights */
              e_st.site = &r_emlrtRSI;
              f_st.site = &td_emlrtRSI;
              b_dynamic_size_checks(&f_st, b_Spower, wt, b_Spower->size[1],
                                    wt->size[0]);
              i1 = b_Spower_iter->size[0];
              b_Spower_iter->size[0] = b_Spower->size[0];
              emxEnsureCapacity_real32_T(&e_st, b_Spower_iter, i1, &yb_emlrtRTEI);
              loop_ub = b_Spower->size[0];
              for (i1 = 0; i1 < loop_ub; i1++) {
                b_Spower_iter->data[i1] = 0.0F;
                b_loop_ub = b_Spower->size[1];
                for (i2 = 0; i2 < b_loop_ub; i2++) {
                  b_Spower_iter->data[i1] += b_Spower->data[i1 + b_Spower->size
                    [0] * i2] * (real32_T)wt->data[i2];
                }
              }
            }

            /* Append the spectrum to the spectrogram */
            i1 = n + 1;
            if ((i1 < 1) || (i1 > mt_spectrogram->size[1])) {
              emlrtDynamicBoundsCheckR2012b(i1, 1, mt_spectrogram->size[1],
                &b_emlrtBCI, &d_st);
            }

            b_loop_ub = freq_inds->size[1];
            for (c_loop_ub = 0; c_loop_ub < b_loop_ub; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub] && ((c_loop_ub + 1 < 1) ||
                   (c_loop_ub + 1 > b_Spower_iter->size[0]))) {
                emlrtDynamicBoundsCheckR2012b(c_loop_ub + 1, 1,
                  b_Spower_iter->size[0], &emlrtBCI, &d_st);
              }
            }

            b_loop_ub = freq_inds->size[1] - 1;
            loop_ub = 0;
            for (c_loop_ub = 0; c_loop_ub <= b_loop_ub; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                loop_ub++;
              }
            }

            i1 = r2->size[0] * r2->size[1];
            r2->size[0] = 1;
            r2->size[1] = loop_ub;
            emxEnsureCapacity_int32_T(&d_st, r2, i1, &db_emlrtRTEI);
            loop_ub = 0;
            for (c_loop_ub = 0; c_loop_ub <= b_loop_ub; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                r2->data[loop_ub] = c_loop_ub + 1;
                loop_ub++;
              }
            }

            emlrtSubAssignSizeCheckR2012b(&mt_spectrogram->size[0], 1, &r2->
              size[1], 1, &emlrtECI, &d_st);
            b_loop_ub = freq_inds->size[1] - 1;
            loop_ub = 0;
            for (c_loop_ub = 0; c_loop_ub <= b_loop_ub; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                loop_ub++;
              }
            }

            i1 = r1->size[0] * r1->size[1];
            r1->size[0] = 1;
            r1->size[1] = loop_ub;
            emxEnsureCapacity_int32_T(&d_st, r1, i1, &db_emlrtRTEI);
            loop_ub = 0;
            for (c_loop_ub = 0; c_loop_ub <= b_loop_ub; c_loop_ub++) {
              if (freq_inds->data[c_loop_ub]) {
                r1->data[loop_ub] = c_loop_ub + 1;
                loop_ub++;
              }
            }

            loop_ub = r1->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              mt_spectrogram->data[i1 + mt_spectrogram->size[0] * n] =
                b_Spower_iter->data[r1->data[i1] - 1];
            }
          }
        }

        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(&d_st);
        }
      } else {
        emlrtHadParallelError = true;
      }
    }

    if (!emlrtHadParallelError) {
      emlrtHeapReferenceStackLeaveScope(&d_st, 15);
      emxFree_real32_T(&wk);
      emxFree_real32_T(&b_b);
      emxFree_real32_T(&b_Spower_iter);
      emxFree_real32_T(&a);
      emxFree_real32_T(&b_Spower);
      emxFree_creal32_T(&fft_data);
      emxFree_real_T(&r5);
      emxFree_boolean_T(&Spower_iter);
      emxFree_real32_T(&Spower);
      emxFree_real32_T(&b_y);
      emxFree_real_T(&r4);
      emxFree_real32_T(&r3);
      emxFree_real32_T(&y);
      emxFree_int32_T(&r2);
      emxFree_int32_T(&r1);
    }
  }

  emlrtPopJmpBuf(sp, &emlrtJBStack);
  emlrtExitParallelRegion(sp, omp_in_parallel());
  emxFree_real_T(&wt);
  emxFree_boolean_T(&freq_inds);

  /* Compute one-sided PSD spectrum  */
  st.site = &s_emlrtRSI;
  i = b_sfreqs->size[0] * b_sfreqs->size[1];
  b_sfreqs->size[0] = 1;
  b_sfreqs->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(&st, b_sfreqs, i, &jb_emlrtRTEI);
  input_sizes_idx_0 = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    b_sfreqs->data[i] = (sfreqs->data[i] == 0.0);
  }

  emxInit_int32_T(&st, &ii, 2, &sc_emlrtRTEI, true);
  b_st.site = &gf_emlrtRSI;
  eml_find(&b_st, b_sfreqs, ii);
  i = window_start->size[0] * window_start->size[1];
  window_start->size[0] = 1;
  window_start->size[1] = ii->size[1];
  emxEnsureCapacity_real_T(&st, window_start, i, &kb_emlrtRTEI);
  input_sizes_idx_0 = ii->size[0] * ii->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    window_start->data[i] = ii->data[i];
  }

  st.site = &t_emlrtRSI;
  b = Fs / 2.0;
  i = b_sfreqs->size[0] * b_sfreqs->size[1];
  b_sfreqs->size[0] = 1;
  b_sfreqs->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(&st, b_sfreqs, i, &lb_emlrtRTEI);
  input_sizes_idx_0 = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    b_sfreqs->data[i] = (sfreqs->data[i] == b);
  }

  emxInit_real_T(&st, &c_y, 2, &tc_emlrtRTEI, true);
  b_st.site = &gf_emlrtRSI;
  eml_find(&b_st, b_sfreqs, ii);
  emxFree_boolean_T(&b_sfreqs);
  if (sfreqs->size[1] < 1) {
    c_y->size[0] = 1;
    c_y->size[1] = 0;
  } else {
    i = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = (int32_T)((real_T)sfreqs->size[1] - 1.0) + 1;
    emxEnsureCapacity_real_T(sp, c_y, i, &nb_emlrtRTEI);
    input_sizes_idx_0 = (int32_T)((real_T)sfreqs->size[1] - 1.0);
    for (i = 0; i <= input_sizes_idx_0; i++) {
      c_y->data[i] = (real_T)i + 1.0;
    }
  }

  emxInit_real_T(sp, &b_window_start, 2, &pb_emlrtRTEI, true);
  st.site = &u_emlrtRSI;
  b_st.site = &jf_emlrtRSI;
  i = b_window_start->size[0] * b_window_start->size[1];
  b_window_start->size[0] = 1;
  b_window_start->size[1] = window_start->size[1] + ii->size[1];
  emxEnsureCapacity_real_T(&b_st, b_window_start, i, &pb_emlrtRTEI);
  input_sizes_idx_0 = window_start->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    b_window_start->data[i] = window_start->data[i];
  }

  input_sizes_idx_0 = ii->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    b_window_start->data[i + window_start->size[1]] = ii->data[i];
  }

  emxInit_real_T(&b_st, &b_select, 2, &rc_emlrtRTEI, true);
  emxInit_real32_T(&b_st, &c_b, 2, &sb_emlrtRTEI, true);
  emxInit_int32_T(&b_st, &ia, 1, &db_emlrtRTEI, true);
  c_st.site = &kf_emlrtRSI;
  do_vectors(&c_st, c_y, b_window_start, b_select, ia, ib_size);
  input_sizes_idx_0 = mt_spectrogram->size[1];
  i = c_b->size[0] * c_b->size[1];
  c_b->size[0] = b_select->size[1];
  c_b->size[1] = mt_spectrogram->size[1];
  emxEnsureCapacity_real32_T(sp, c_b, i, &sb_emlrtRTEI);
  emxFree_real_T(&b_window_start);
  emxFree_int32_T(&ia);
  emxFree_real_T(&c_y);
  for (i = 0; i < input_sizes_idx_0; i++) {
    nz = b_select->size[1];
    for (b_i = 0; b_i < nz; b_i++) {
      if (b_select->data[b_i] != (int32_T)muDoubleScalarFloor(b_select->data[b_i]))
      {
        emlrtIntegerCheckR2012b(b_select->data[b_i], &e_emlrtDCI, sp);
      }

      b_input_sizes_idx_0 = (int32_T)b_select->data[b_i];
      if ((b_input_sizes_idx_0 < 1) || (b_input_sizes_idx_0 >
           mt_spectrogram->size[0])) {
        emlrtDynamicBoundsCheckR2012b(b_input_sizes_idx_0, 1,
          mt_spectrogram->size[0], &g_emlrtBCI, sp);
      }

      c_b->data[b_i + c_b->size[0] * i] = mt_spectrogram->data
        [(b_input_sizes_idx_0 + mt_spectrogram->size[0] * i) - 1];
    }
  }

  emxFree_real_T(&b_select);
  input_sizes_idx_0 = c_b->size[0] * c_b->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    c_b->data[i] *= 2.0F;
  }

  emxInit_real32_T(sp, &varargin_1, 2, &ub_emlrtRTEI, true);
  st.site = &v_emlrtRSI;
  input_sizes_idx_0 = mt_spectrogram->size[1];
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = window_start->size[1];
  varargin_1->size[1] = mt_spectrogram->size[1];
  emxEnsureCapacity_real32_T(&st, varargin_1, i, &ub_emlrtRTEI);
  for (i = 0; i < input_sizes_idx_0; i++) {
    nz = window_start->size[1];
    for (b_i = 0; b_i < nz; b_i++) {
      b_input_sizes_idx_0 = (int32_T)window_start->data[b_i];
      if ((b_input_sizes_idx_0 < 1) || (b_input_sizes_idx_0 >
           mt_spectrogram->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)window_start->data[b_i], 1,
          mt_spectrogram->size[0], &h_emlrtBCI, &st);
      }

      varargin_1->data[b_i + varargin_1->size[0] * i] = mt_spectrogram->data
        [(b_input_sizes_idx_0 + mt_spectrogram->size[0] * i) - 1];
    }
  }

  emxInit_real32_T(&st, &varargin_3, 2, &vb_emlrtRTEI, true);
  input_sizes_idx_0 = mt_spectrogram->size[1];
  i = varargin_3->size[0] * varargin_3->size[1];
  varargin_3->size[0] = ii->size[1];
  varargin_3->size[1] = mt_spectrogram->size[1];
  emxEnsureCapacity_real32_T(&st, varargin_3, i, &vb_emlrtRTEI);
  for (i = 0; i < input_sizes_idx_0; i++) {
    nz = ii->size[1];
    for (b_i = 0; b_i < nz; b_i++) {
      if ((ii->data[b_i] < 1) || (ii->data[b_i] > mt_spectrogram->size[0])) {
        emlrtDynamicBoundsCheckR2012b(ii->data[b_i], 1, mt_spectrogram->size[0],
          &i_emlrtBCI, &st);
      }

      varargin_3->data[b_i + varargin_3->size[0] * i] = mt_spectrogram->data
        [(ii->data[b_i] + mt_spectrogram->size[0] * i) - 1];
    }
  }

  b_st.site = &wf_emlrtRSI;
  if ((window_start->size[1] != 0) && (mt_spectrogram->size[1] != 0)) {
    nz = mt_spectrogram->size[1];
  } else if ((c_b->size[0] != 0) && (c_b->size[1] != 0)) {
    nz = c_b->size[1];
  } else if ((ii->size[1] != 0) && (mt_spectrogram->size[1] != 0)) {
    nz = mt_spectrogram->size[1];
  } else {
    if (mt_spectrogram->size[1] > 0) {
      nz = mt_spectrogram->size[1];
    } else {
      nz = 0;
    }

    if (c_b->size[1] > nz) {
      nz = c_b->size[1];
    }

    if (mt_spectrogram->size[1] > nz) {
      nz = mt_spectrogram->size[1];
    }
  }

  c_st.site = &xf_emlrtRSI;
  if ((mt_spectrogram->size[1] != nz) && ((window_start->size[1] != 0) &&
       (mt_spectrogram->size[1] != 0))) {
    emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
      "MATLAB:catenate:matrixDimensionMismatch",
      "MATLAB:catenate:matrixDimensionMismatch", 0);
  }

  if ((c_b->size[1] != nz) && ((c_b->size[0] != 0) && (c_b->size[1] != 0))) {
    emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
      "MATLAB:catenate:matrixDimensionMismatch",
      "MATLAB:catenate:matrixDimensionMismatch", 0);
  }

  if ((mt_spectrogram->size[1] != nz) && ((ii->size[1] != 0) &&
       (mt_spectrogram->size[1] != 0))) {
    emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
      "MATLAB:catenate:matrixDimensionMismatch",
      "MATLAB:catenate:matrixDimensionMismatch", 0);
  }

  empty_non_axis_sizes = (nz == 0);
  if (empty_non_axis_sizes || ((window_start->size[1] != 0) &&
       (mt_spectrogram->size[1] != 0))) {
    b_input_sizes_idx_0 = window_start->size[1];
  } else {
    b_input_sizes_idx_0 = 0;
  }

  emxFree_real_T(&window_start);
  if (empty_non_axis_sizes || ((c_b->size[0] != 0) && (c_b->size[1] != 0))) {
    input_sizes_idx_0 = c_b->size[0];
  } else {
    input_sizes_idx_0 = 0;
  }

  if (empty_non_axis_sizes || ((ii->size[1] != 0) && (mt_spectrogram->size[1] !=
        0))) {
    winsize_samples = ii->size[1];
  } else {
    winsize_samples = 0;
  }

  emxFree_int32_T(&ii);
  num_tapers = b_input_sizes_idx_0;
  i = mt_spectrogram->size[0] * mt_spectrogram->size[1];
  mt_spectrogram->size[0] = (num_tapers + input_sizes_idx_0) + winsize_samples;
  mt_spectrogram->size[1] = nz;
  emxEnsureCapacity_real32_T(sp, mt_spectrogram, i, &cc_emlrtRTEI);
  for (i = 0; i < nz; i++) {
    for (b_i = 0; b_i < num_tapers; b_i++) {
      mt_spectrogram->data[b_i + mt_spectrogram->size[0] * i] = varargin_1->
        data[b_i + num_tapers * i] / (real32_T)Fs;
    }
  }

  emxFree_real32_T(&varargin_1);
  for (i = 0; i < nz; i++) {
    for (b_i = 0; b_i < input_sizes_idx_0; b_i++) {
      mt_spectrogram->data[(b_i + num_tapers) + mt_spectrogram->size[0] * i] =
        c_b->data[b_i + input_sizes_idx_0 * i] / (real32_T)Fs;
    }
  }

  emxFree_real32_T(&c_b);
  for (i = 0; i < nz; i++) {
    for (b_i = 0; b_i < winsize_samples; b_i++) {
      mt_spectrogram->data[((b_i + num_tapers) + input_sizes_idx_0) +
        mt_spectrogram->size[0] * i] = varargin_3->data[b_i + winsize_samples *
        i] / (real32_T)Fs;
    }
  }

  emxFree_real32_T(&varargin_3);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (multitaper_spectrogram_coder.c) */
