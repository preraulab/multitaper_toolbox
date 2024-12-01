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
#include "assertCompatibleDims.h"
#include "colon.h"
#include "combineVectorElements.h"
#include "detrend.h"
#include "div.h"
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
#include "repmat.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "mwmathutil.h"
#include "omp.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    55,                             /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    61,                             /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    65,                             /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    77,                             /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    98,                             /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    101,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    105,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    112,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    114,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo j_emlrtRSI = {
    118,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI = {
    121,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo l_emlrtRSI = {
    124,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo m_emlrtRSI = {
    129,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI = {
    130,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo o_emlrtRSI = {
    134,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo p_emlrtRSI = {
    136,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo q_emlrtRSI = {
    137,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI = {
    142,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI = {
    150,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo t_emlrtRSI = {
    151,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo u_emlrtRSI = {
    152,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo v_emlrtRSI = {
    153,                            /* lineNo */
    "multitaper_spectrogram_coder", /* fcnName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pathName */
};

static emlrtRSInfo w_emlrtRSI = {
    125,     /* lineNo */
    "colon", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo db_emlrtRSI = {
    15,    /* lineNo */
    "sum", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/datafun/sum.m" /* pathName
                                                                            */
};

static emlrtRSInfo nd_emlrtRSI = {
    69,                  /* lineNo */
    "eml_mtimes_helper", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/"
    "eml_mtimes_helper.m" /* pathName */
};

static emlrtRSInfo td_emlrtRSI = {
    34,               /* lineNo */
    "rdivide_helper", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
    "rdivide_helper.m" /* pathName */
};

static emlrtRSInfo
    ud_emlrtRSI =
        {
            53,    /* lineNo */
            "div", /* fcnName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
            "div.m" /* pathName */
};

static emlrtRSInfo oe_emlrtRSI = {
    39,     /* lineNo */
    "find", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                           */
};

static emlrtRSInfo re_emlrtRSI = {
    19,        /* lineNo */
    "setdiff", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/setdiff.m" /* pathName
                                                                            */
};

static emlrtRSInfo se_emlrtRSI = {
    97,          /* lineNo */
    "eml_setop", /* fcnName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/private/"
    "eml_setop.m" /* pathName */
};

static emlrtRSInfo
    xe_emlrtRSI =
        {
            41,    /* lineNo */
            "cat", /* fcnName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
            "cat.m" /* pathName */
};

static emlrtRSInfo
    ye_emlrtRSI =
        {
            65,         /* lineNo */
            "cat_impl", /* fcnName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
            "cat.m" /* pathName */
};

static emlrtRTEInfo
    b_emlrtRTEI =
        {
            225,                   /* lineNo */
            27,                    /* colNo */
            "check_non_axis_size", /* fName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
            "cat.m" /* pName */
};

static emlrtECInfo emlrtECI = {
    -1,                             /* nDims */
    146,                            /* lineNo */
    5,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtBCInfo emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    146,                            /* lineNo */
    22,                             /* colNo */
    "mt_spectrogram",               /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    106,                            /* lineNo */
    26,                             /* colNo */
    "mt_spectrogram",               /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtECInfo b_emlrtECI = {
    2,                              /* nDims */
    137,                            /* lineNo */
    29,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtECInfo c_emlrtECI = {
    1,                              /* nDims */
    137,                            /* lineNo */
    29,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtECInfo d_emlrtECI = {
    2,                              /* nDims */
    136,                            /* lineNo */
    16,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtECInfo e_emlrtECI = {
    1,                              /* nDims */
    136,                            /* lineNo */
    16,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtDCInfo emlrtDCI = {
    136,                            /* lineNo */
    30,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    1              /* checkKind */
};

static emlrtECInfo f_emlrtECI = {
    1,                              /* nDims */
    134,                            /* lineNo */
    50,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtDCInfo b_emlrtDCI = {
    134,                            /* lineNo */
    79,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    1              /* checkKind */
};

static emlrtECInfo g_emlrtECI = {
    2,                              /* nDims */
    124,                            /* lineNo */
    14,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtECInfo h_emlrtECI = {
    2,                              /* nDims */
    118,                            /* lineNo */
    20,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtECInfo i_emlrtECI = {
    1,                              /* nDims */
    118,                            /* lineNo */
    20,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtDCInfo c_emlrtDCI = {
    77,                             /* lineNo */
    24,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    4              /* checkKind */
};

static emlrtECInfo j_emlrtECI = {
    2,                              /* nDims */
    68,                             /* lineNo */
    13,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    69,                             /* lineNo */
    17,                             /* colNo */
    "sfreqs",                       /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    98,                             /* lineNo */
    38,                             /* colNo */
    "window_start",                 /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    98,                             /* lineNo */
    25,                             /* colNo */
    "data",                         /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = {
    98,                             /* lineNo */
    25,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    1              /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    153,                            /* lineNo */
    65,                             /* colNo */
    "mt_spectrogram",               /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = {
    153,                            /* lineNo */
    65,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    1              /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    153,                            /* lineNo */
    34,                             /* colNo */
    "mt_spectrogram",               /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    153,                            /* lineNo */
    91,                             /* colNo */
    "mt_spectrogram",               /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    130,                            /* lineNo */
    37,                             /* colNo */
    "Spower",                       /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = {
    -1,                             /* iFirst */
    -1,                             /* iLast */
    146,                            /* lineNo */
    39,                             /* colNo */
    "mt_spectrum",                  /* aName */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m", /* pName */
    0              /* checkKind */
};

static emlrtRTEInfo x_emlrtRTEI = {
    55,                             /* lineNo */
    1,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo y_emlrtRTEI = {
    65,                             /* lineNo */
    1,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = {
    68,                             /* lineNo */
    14,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = {
    68,                             /* lineNo */
    47,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = {
    68,                             /* lineNo */
    1,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo db_emlrtRTEI = {
    1,                              /* lineNo */
    45,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = {
    73,                             /* lineNo */
    1,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = {
    77,                             /* lineNo */
    1,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = {
    92,                             /* lineNo */
    5,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = {
    90,                             /* lineNo */
    5,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = {
    88,                             /* lineNo */
    5,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo jb_emlrtRTEI = {
    150,                            /* lineNo */
    18,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo kb_emlrtRTEI = {
    150,                            /* lineNo */
    1,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = {
    151,                            /* lineNo */
    23,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo mb_emlrtRTEI = {
    98,                             /* lineNo */
    44,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo nb_emlrtRTEI = {
    28,      /* lineNo */
    9,       /* colNo */
    "colon", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/colon.m" /* pName
                                                                          */
};

static emlrtRTEInfo ob_emlrtRTEI = {
    98,                             /* lineNo */
    5,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = {
    152,                            /* lineNo */
    36,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = {
    101,                            /* lineNo */
    12,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = {
    16,      /* lineNo */
    13,      /* colNo */
    "isnan", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/isnan.m" /* pName
                                                                            */
};

static emlrtRTEInfo sb_emlrtRTEI = {
    153,                            /* lineNo */
    50,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo tb_emlrtRTEI = {
    112,                            /* lineNo */
    32,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo ub_emlrtRTEI = {
    153,                            /* lineNo */
    19,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = {
    153,                            /* lineNo */
    76,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = {
    124,                            /* lineNo */
    14,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = {
    124,                            /* lineNo */
    34,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo yb_emlrtRTEI = {
    142,                            /* lineNo */
    9,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo ac_emlrtRTEI = {
    130,                            /* lineNo */
    28,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo bc_emlrtRTEI = {
    153,                            /* lineNo */
    1,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo cc_emlrtRTEI = {
    131,                            /* lineNo */
    9,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo dc_emlrtRTEI = {
    146,                            /* lineNo */
    27,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo ec_emlrtRTEI = {
    134,                            /* lineNo */
    50,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo fc_emlrtRTEI = {
    134,                            /* lineNo */
    74,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo gc_emlrtRTEI = {
    88,                  /* lineNo */
    9,                   /* colNo */
    "eml_mtimes_helper", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/ops/"
    "eml_mtimes_helper.m" /* pName */
};

static emlrtRTEInfo hc_emlrtRTEI = {
    136,                            /* lineNo */
    25,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo ic_emlrtRTEI = {
    137,                            /* lineNo */
    29,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo jc_emlrtRTEI = {
    137,                            /* lineNo */
    25,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo kc_emlrtRTEI = {
    136,                            /* lineNo */
    13,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo lc_emlrtRTEI = {
    130,                            /* lineNo */
    9,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo mc_emlrtRTEI = {
    124,                            /* lineNo */
    5,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo nc_emlrtRTEI = {
    121,                            /* lineNo */
    5,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo oc_emlrtRTEI = {
    98,                             /* lineNo */
    20,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo pc_emlrtRTEI = {
    136,                            /* lineNo */
    16,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo qc_emlrtRTEI = {
    137,                            /* lineNo */
    47,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo rc_emlrtRTEI = {
    146,                            /* lineNo */
    39,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo sc_emlrtRTEI = {
    152,                            /* lineNo */
    1,                              /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo tc_emlrtRTEI = {
    68,                             /* lineNo */
    13,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo uc_emlrtRTEI = {
    31,     /* lineNo */
    6,      /* colNo */
    "find", /* fName */
    "/Applications/MATLAB_R2024b.app/toolbox/eml/lib/matlab/elmat/find.m" /* pName
                                                                           */
};

static emlrtRTEInfo vc_emlrtRTEI = {
    152,                            /* lineNo */
    18,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRTEInfo ee_emlrtRTEI = {
    118,                            /* lineNo */
    20,                             /* colNo */
    "multitaper_spectrogram_coder", /* fName */
    "/Users/alexhe/Dropbox "
    "(Personal)/Active_projects/preraulab_code/multitaper/matlab/"
    "multitaper_spectrogram_mex/multitaper_spectrog"
    "ram_coder.m" /* pName */
};

static emlrtRSInfo
    bf_emlrtRSI =
        {
            54,    /* lineNo */
            "div", /* fcnName */
            "/Applications/MATLAB_R2024b.app/toolbox/eml/eml/+coder/+internal/"
            "div.m" /* pathName */
};

/* Function Declarations */
static void binary_expand_op(const emlrtStack *sp, emxArray_real32_T *in1,
                             const emlrtRSInfo in2,
                             const emxArray_real32_T *in3,
                             const emxArray_real32_T *in4);

static void binary_expand_op_1(const emlrtStack *sp, emxArray_real32_T *in1,
                               const emxArray_real_T *in2);

static void binary_expand_op_2(const emlrtStack *sp, emxArray_real32_T *in1,
                               const emxArray_real_T *in2);

static void plus(const emlrtStack *sp, emxArray_real32_T *in1,
                 const emxArray_real32_T *in2);

/* Function Definitions */
static void binary_expand_op(const emlrtStack *sp, emxArray_real32_T *in1,
                             const emlrtRSInfo in2,
                             const emxArray_real32_T *in3,
                             const emxArray_real32_T *in4)
{
  emlrtStack st;
  emxArray_real32_T *b_in3;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T b_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_0_1;
  int32_T stride_1_0;
  int32_T stride_1_1;
  const real32_T *in3_data;
  const real32_T *in4_data;
  real32_T *b_in3_data;
  st.prev = sp;
  st.tls = sp->tls;
  in4_data = in4->data;
  in3_data = in3->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real32_T(sp, &b_in3, 2, &ic_emlrtRTEI);
  if (in4->size[1] == 1) {
    loop_ub = in3->size[1];
  } else {
    loop_ub = in4->size[1];
  }
  i = b_in3->size[0] * b_in3->size[1];
  b_in3->size[0] = loop_ub;
  if (in4->size[0] == 1) {
    b_loop_ub = in3->size[0];
  } else {
    b_loop_ub = in4->size[0];
  }
  b_in3->size[1] = b_loop_ub;
  emxEnsureCapacity_real32_T(sp, b_in3, i, &ic_emlrtRTEI);
  b_in3_data = b_in3->data;
  stride_0_0 = (in3->size[1] != 1);
  stride_0_1 = (in3->size[0] != 1);
  stride_1_0 = (in4->size[1] != 1);
  stride_1_1 = (in4->size[0] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in3_data[i1 + b_in3->size[0] * i] =
          in3_data[aux_0_1 + in3->size[0] * (i1 * stride_0_0)] *
          in4_data[aux_1_1 + in4->size[0] * (i1 * stride_1_0)];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  st.site = (emlrtRSInfo *)&in2;
  sum(&st, b_in3, in1);
  emxFree_real32_T(sp, &b_in3);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_1(const emlrtStack *sp, emxArray_real32_T *in1,
                               const emxArray_real_T *in2)
{
  emxArray_real32_T *b_in1;
  const real_T *in2_data;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T b_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_0_1;
  int32_T stride_1_0;
  int32_T stride_1_1;
  real32_T *b_in1_data;
  real32_T *in1_data;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real32_T(sp, &b_in1, 2, &pc_emlrtRTEI);
  if (in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = in2->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  if (in2->size[1] == 1) {
    b_loop_ub = in1->size[1];
  } else {
    b_loop_ub = in2->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real32_T(sp, b_in1, i, &pc_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_0_1 = (in1->size[1] != 1);
  stride_1_0 = (in2->size[0] != 1);
  stride_1_1 = (in2->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * aux_0_1] *
          (real32_T)in2_data[i1 * stride_1_0 + in2->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = loop_ub;
  in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real32_T(sp, in1, i, &pc_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real32_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_2(const emlrtStack *sp, emxArray_real32_T *in1,
                               const emxArray_real_T *in2)
{
  emxArray_real32_T *b_in1;
  const real_T *in2_data;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T b_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_0_1;
  int32_T stride_1_0;
  int32_T stride_1_1;
  real32_T *b_in1_data;
  real32_T *in1_data;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real32_T(sp, &b_in1, 2, &ee_emlrtRTEI);
  if (in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = in2->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  if (in2->size[1] == 1) {
    b_loop_ub = in1->size[1];
  } else {
    b_loop_ub = in2->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real32_T(sp, b_in1, i, &ee_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_0_1 = (in1->size[1] != 1);
  stride_1_0 = (in2->size[0] != 1);
  stride_1_1 = (in2->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * aux_0_1] *
          (real32_T)in2_data[i1 * stride_1_0 + in2->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = loop_ub;
  in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real32_T(sp, in1, i, &ee_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real32_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void plus(const emlrtStack *sp, emxArray_real32_T *in1,
                 const emxArray_real32_T *in2)
{
  emxArray_real32_T *b_in1;
  int32_T b_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_1_0;
  const real32_T *in2_data;
  real32_T *b_in1_data;
  real32_T *in1_data;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real32_T(sp, &b_in1, 2, &ec_emlrtRTEI);
  if (in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = in2->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  b_loop_ub = in1->size[1];
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real32_T(sp, b_in1, i, &ec_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_1_0 = (in2->size[0] != 1);
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * i] +
          in2_data[i1 * stride_1_0 + in2->size[0] * i];
    }
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = loop_ub;
  emxEnsureCapacity_real32_T(sp, in1, i, &ec_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real32_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

emlrtCTX emlrtGetRootTLSGlobal(void)
{
  return emlrtRootTLSGlobal;
}

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, emlrtConstCTX aTLS,
                         void *aData)
{
  omp_set_lock(&emlrtLockGlobal);
  emlrtCallLockeeFunction(aLockee, aTLS, aData);
  omp_unset_lock(&emlrtLockGlobal);
}

void multitaper_spectrogram_coder(
    const emlrtStack *sp, const emxArray_real32_T *data, real_T Fs,
    const real_T frequency_range[2], const emxArray_real_T *DPSS_tapers,
    const emxArray_real_T *DPSS_eigen, real_T winstep_samples, real_T min_NFFT,
    real_T detrend_opt, real_T weighting, emxArray_real32_T *mt_spectrogram,
    emxArray_real_T *stimes, emxArray_real_T *sfreqs)
{
  jmp_buf emlrtJBEnviron;
  jmp_buf *volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  emxArray_boolean_T *a;
  emxArray_boolean_T *freq_inds;
  emxArray_boolean_T *r;
  emxArray_creal32_T *fft_data;
  emxArray_int32_T *ia;
  emxArray_int32_T *ii;
  emxArray_int32_T *r5;
  emxArray_int32_T *r8;
  emxArray_real32_T *Spower;
  emxArray_real32_T *Spower_iter;
  emxArray_real32_T *b_Spower;
  emxArray_real32_T *b_a;
  emxArray_real32_T *b_b;
  emxArray_real32_T *b_wk;
  emxArray_real32_T *b_y;
  emxArray_real32_T *r6;
  emxArray_real32_T *varargin_1;
  emxArray_real32_T *varargin_3;
  emxArray_real32_T *wk;
  emxArray_real32_T *y;
  emxArray_real_T *b_select;
  emxArray_real_T *b_window_start;
  emxArray_real_T *c_y;
  emxArray_real_T *r7;
  emxArray_real_T *window_start;
  emxArray_real_T *wt;
  creal32_T *fft_data_data;
  const real_T *DPSS_eigen_data;
  const real_T *DPSS_tapers_data;
  real_T b;
  real_T c_window_start;
  real_T d;
  real_T nfft;
  real_T *r3;
  real_T *sfreqs_data;
  real_T *stimes_data;
  real_T *window_start_data;
  int32_T b_i;
  int32_T b_ii;
  int32_T b_input_sizes_idx_0;
  int32_T b_loop_ub;
  int32_T c_loop_ub;
  int32_T d_loop_ub;
  int32_T e_loop_ub;
  int32_T f_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T input_sizes_idx_0;
  int32_T loop_ub;
  int32_T multitaper_spectrogram_coder_numThreads;
  int32_T n;
  int32_T num_tapers;
  int32_T num_tapers_tmp;
  int32_T nz_tmp;
  int32_T sizes_idx_0;
  int32_T winsize_samples;
  int32_T *ii_data;
  int32_T *r2;
  int32_T *r4;
  const real32_T *data_data;
  real32_T Tpower;
  real32_T *Spower_data;
  real32_T *Spower_iter_data;
  real32_T *a_data;
  real32_T *b_data;
  real32_T *b_wk_data;
  real32_T *mt_spectrogram_data;
  real32_T *varargin_1_data;
  real32_T *varargin_3_data;
  real32_T *wk_data;
  real32_T *y_data;
  boolean_T c_b;
  boolean_T emlrtHadParallelError = false;
  boolean_T empty_non_axis_sizes;
  boolean_T *b_a_data;
  boolean_T *freq_inds_data;
  boolean_T *r1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  sfreqs_data = sfreqs->data;
  DPSS_eigen_data = DPSS_eigen->data;
  DPSS_tapers_data = DPSS_tapers->data;
  data_data = data->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /* MULTITAPER_SPECTROGRAM  Compute the multitaper spectrogram for time series
   * data */
  /*  */
  /*    This is the coder portion for mex compilation. It takes processed */
  /*    multitaper parameters from multitaper_spectrogram_coder_mex.m as inputs
   */
  /*    and skips internal input processing.  */
  /*     */
  /*    Usage: */
  /*    Direct input: */
  /*    [spect,stimes,sfreqs] = multitaper_spectrogram_coder(data, Fs,
   * frequency_range, DPSS_tapers, DPSS_eigen, winstep_samples, min_NFFT,
   * detrend_opt, weighting) */
  /*  */
  /*    Input: */
  /*    data: <number of samples> x 1 vector - time series data -- required */
  /*    Fs: double - sampling frequency in Hz  -- required */
  /*    frequency_range: 1x2 vector - [<min frequency>, <max frequency>] --
   * required  */
  /*    DPSS_tapers: Nxk matrix - Slepian tapers -- required  */
  /*    DPSS_eigen: 1xk vector - eigenvalues of the Slepian tapers -- required
   */
  /*    winstep_samples: double - number of samples in step size of windows --
   * required  */
  /*    min_NFFT: double - minimum allowable NFFT size, adds zero padding for
   * interpolation (closest 2^x) -- required */
  /*    detrend_opt: double - how to detrend data window (2='linear',
   * 1='constant', 0='off') -- required */
  /*    weighting: double - how to weight the tapers (0='unity', 1='eigen',
   * 2='adapt') -- required */
  /*  */
  /*    Output: */
  /*    mt_spectrogram: FxT matrix of one-sided power spectral density (PSD) */
  /*    stimes: 1XT vector of times for the center of the spectral bins */
  /*    sfreqs: 1XF vector of frequency bins for the spectrogram */
  /*  */
  /*    This code is companion to the paper: */
  /*          "Sleep Neurophysiological Dynamics Through the Lens of Multitaper
   * Spectral Analysis" */
  /*          Michael J. Prerau, Ritchie E. Brown, Matt T. Bianchi, Jeffrey M.
   * Ellenbogen, Patrick L. Purdon */
  /*          December 7, 2016 : 60-92 */
  /*          DOI: 10.1152/physiol.00062.2015 */
  /*    which should be cited for academic use of this code. */
  /*  */
  /*    A full tutorial on the multitaper spectrogram can be found at: */
  /*    http://www.sleepEEG.org/multitaper */
  /*  */
  /*    Copyright 2024 Michael J. Prerau, Ph.D. - http://www.sleepEEG.org */
  /*    This work is licensed under a Creative Commons
   * Attribution-NonCommercial-ShareAlike 4.0 International License. */
  /*    (http://creativecommons.org/licenses/by-nc-sa/4.0/) */
  /*  */
  /*    Last modified 2/16/2021 */
  /*  ******************************************************************** */
  /* Generate DPSS tapers (STEP 1) */
  /*  Done outside this _coder function. */
  /* Get taper matrix dimensions  */
  winsize_samples = DPSS_tapers->size[0];
  num_tapers_tmp = DPSS_tapers->size[1];
  num_tapers = DPSS_tapers->size[1];
  /* Total data length */
  /* Window start indices */
  st.site = &emlrtRSI;
  b = (real_T)(data->size[0] - DPSS_tapers->size[0]) + 1.0;
  emxInit_real_T(&st, &window_start, 2, &x_emlrtRTEI);
  window_start_data = window_start->data;
  if (muDoubleScalarIsNaN(winstep_samples)) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &x_emlrtRTEI);
    window_start_data = window_start->data;
    window_start_data[0] = rtNaN;
  } else if ((winstep_samples == 0.0) ||
             ((b > 1.0) && (winstep_samples < 0.0)) ||
             ((b < 1.0) && (winstep_samples > 0.0))) {
    window_start->size[0] = 1;
    window_start->size[1] = 0;
  } else if (muDoubleScalarIsInf(winstep_samples)) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &x_emlrtRTEI);
    window_start_data = window_start->data;
    window_start_data[0] = 1.0;
  } else if (muDoubleScalarFloor(winstep_samples) == winstep_samples) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    loop_ub = (int32_T)((b - 1.0) / winstep_samples);
    window_start->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &x_emlrtRTEI);
    window_start_data = window_start->data;
    for (i = 0; i <= loop_ub; i++) {
      window_start_data[i] = winstep_samples * (real_T)i + 1.0;
    }
  } else {
    b_st.site = &w_emlrtRSI;
    eml_float_colon(&b_st, winstep_samples, b, window_start);
    window_start_data = window_start->data;
  }
  /* Number of windows */
  /* Number of points in the FFT */
  st.site = &b_emlrtRSI;
  b_st.site = &bb_emlrtRSI;
  st.site = &b_emlrtRSI;
  b_st.site = &bb_emlrtRSI;
  nfft = muDoubleScalarMax(
      muDoubleScalarMax(
          muDoubleScalarPower(2.0, nextpow2(DPSS_tapers->size[0])),
          DPSS_tapers->size[0]),
      muDoubleScalarPower(2.0, nextpow2(min_NFFT)));
  /* Create the frequency vector */
  b = Fs / nfft;
  st.site = &c_emlrtRSI;
  if (muDoubleScalarIsNaN(b) || muDoubleScalarIsNaN(Fs)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &y_emlrtRTEI);
    sfreqs_data = sfreqs->data;
    sfreqs_data[0] = rtNaN;
  } else if ((b == 0.0) || ((Fs > 0.0) && (b < 0.0)) ||
             ((Fs < 0.0) && (b > 0.0))) {
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 0;
  } else if (muDoubleScalarIsInf(Fs) && muDoubleScalarIsInf(b)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &y_emlrtRTEI);
    sfreqs_data = sfreqs->data;
    sfreqs_data[0] = rtNaN;
  } else if (muDoubleScalarIsInf(b)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &y_emlrtRTEI);
    sfreqs_data = sfreqs->data;
    sfreqs_data[0] = 0.0;
  } else if (muDoubleScalarFloor(b) == b) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    loop_ub = (int32_T)(Fs / b);
    sfreqs->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &y_emlrtRTEI);
    sfreqs_data = sfreqs->data;
    for (i = 0; i <= loop_ub; i++) {
      sfreqs_data[i] = b * (real_T)i;
    }
  } else {
    b_st.site = &w_emlrtRSI;
    b_eml_float_colon(&b_st, b, Fs, sfreqs);
    sfreqs_data = sfreqs->data;
  }
  /*  all possible frequencies */
  /* Get just the frequencies for the given frequency range */
  emxInit_boolean_T(sp, &freq_inds, 2, &cb_emlrtRTEI);
  i = freq_inds->size[0] * freq_inds->size[1];
  freq_inds->size[0] = 1;
  loop_ub = sfreqs->size[1];
  freq_inds->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, freq_inds, i, &ab_emlrtRTEI);
  freq_inds_data = freq_inds->data;
  emxInit_boolean_T(sp, &r, 2, &tc_emlrtRTEI);
  i = r->size[0] * r->size[1];
  r->size[0] = 1;
  r->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, r, i, &bb_emlrtRTEI);
  r1 = r->data;
  for (i = 0; i < loop_ub; i++) {
    b = sfreqs_data[i];
    freq_inds_data[i] = (b >= frequency_range[0]);
    r1[i] = (b <= frequency_range[1]);
  }
  if (freq_inds->size[1] != r->size[1]) {
    emlrtSizeEqCheckNDErrorR2021b(&freq_inds->size[0], &r->size[0], &j_emlrtECI,
                                  (emlrtCTX)sp);
  }
  i = freq_inds->size[0] * freq_inds->size[1];
  freq_inds->size[0] = 1;
  emxEnsureCapacity_boolean_T(sp, freq_inds, i, &cb_emlrtRTEI);
  freq_inds_data = freq_inds->data;
  loop_ub = freq_inds->size[1] - 1;
  for (i = 0; i <= loop_ub; i++) {
    freq_inds_data[i] = (freq_inds_data[i] && r1[i]);
  }
  emxFree_boolean_T(sp, &r);
  input_sizes_idx_0 = freq_inds->size[1];
  sizes_idx_0 = 0;
  for (b_i = 0; b_i < input_sizes_idx_0; b_i++) {
    if (freq_inds_data[b_i]) {
      sizes_idx_0++;
    }
  }
  b_input_sizes_idx_0 = 0;
  for (b_i = 0; b_i < input_sizes_idx_0; b_i++) {
    if (freq_inds_data[b_i]) {
      if (b_i > sfreqs->size[1] - 1) {
        emlrtDynamicBoundsCheckR2012b(b_i, 0, sfreqs->size[1] - 1, &c_emlrtBCI,
                                      (emlrtConstCTX)sp);
      }
      sfreqs_data[b_input_sizes_idx_0] = sfreqs_data[b_i];
      b_input_sizes_idx_0++;
    }
  }
  i = sfreqs->size[0] * sfreqs->size[1];
  sfreqs->size[0] = 1;
  sfreqs->size[1] = sizes_idx_0;
  emxEnsureCapacity_real_T(sp, sfreqs, i, &db_emlrtRTEI);
  sfreqs_data = sfreqs->data;
  /* Compute the times of the middle of each spectrum */
  input_sizes_idx_0 =
      (int32_T)muDoubleScalarRound((real_T)DPSS_tapers->size[0] / 2.0);
  i = stimes->size[0] * stimes->size[1];
  stimes->size[0] = 1;
  loop_ub = window_start->size[1];
  stimes->size[1] = window_start->size[1];
  emxEnsureCapacity_real_T(sp, stimes, i, &eb_emlrtRTEI);
  stimes_data = stimes->data;
  for (i = 0; i < loop_ub; i++) {
    stimes_data[i] =
        ((window_start_data[i] + (real_T)input_sizes_idx_0) - 1.0) / Fs;
  }
  /*  stimes start from 0 */
  /* Preallocate spectrogram and slice data for efficient parallel computing */
  st.site = &d_emlrtRSI;
  b_st.site = &db_emlrtRSI;
  c_st.site = &eb_emlrtRSI;
  nz_tmp = combineVectorElements(&c_st, freq_inds);
  if (nz_tmp < 0) {
    emlrtNonNegativeCheckR2012b(nz_tmp, &c_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = mt_spectrogram->size[0] * mt_spectrogram->size[1];
  mt_spectrogram->size[0] = nz_tmp;
  mt_spectrogram->size[1] = window_start->size[1];
  emxEnsureCapacity_real32_T(sp, mt_spectrogram, i, &fb_emlrtRTEI);
  mt_spectrogram_data = mt_spectrogram->data;
  b_i = nz_tmp * window_start->size[1];
  for (i = 0; i < b_i; i++) {
    mt_spectrogram_data[i] = 0.0F;
  }
  /*  COMPUTE THE MULTITAPER SPECTROGRAM */
  /*  */
  /*      STEP 1: Compute DPSS tapers based on desired spectral properties */
  /*      STEP 2: Multiply the data segment by the DPSS Tapers */
  /*      STEP 3: Compute the spectrum for each tapered segment */
  /*      STEP 4: Take the mean of the tapered spectra */
  /*  pre-compute weights */
  emxInit_real_T(sp, &wt, 1, &ib_emlrtRTEI);
  if (weighting == 1.0) {
    sizes_idx_0 = DPSS_eigen->size[0];
    i = wt->size[0];
    wt->size[0] = DPSS_eigen->size[0];
    emxEnsureCapacity_real_T(sp, wt, i, &ib_emlrtRTEI);
    stimes_data = wt->data;
    for (i = 0; i < sizes_idx_0; i++) {
      stimes_data[i] = DPSS_eigen_data[i] / (real_T)DPSS_tapers->size[1];
    }
  } else if (weighting == 0.0) {
    i = wt->size[0];
    wt->size[0] = DPSS_tapers->size[1];
    emxEnsureCapacity_real_T(sp, wt, i, &hb_emlrtRTEI);
    stimes_data = wt->data;
    for (i = 0; i < num_tapers_tmp; i++) {
      stimes_data[i] = 1.0 / (real_T)DPSS_tapers->size[1];
    }
  } else {
    i = wt->size[0];
    wt->size[0] = 1;
    emxEnsureCapacity_real_T(sp, wt, i, &gb_emlrtRTEI);
    stimes_data = wt->data;
    stimes_data[0] = 0.0;
  }
  /* Loop in parallel over all of the windows */
  input_sizes_idx_0 = window_start->size[1];
  emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
  emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
  multitaper_spectrogram_coder_numThreads = emlrtAllocRegionTLSs(
      sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel                                                           \
num_threads(multitaper_spectrogram_coder_numThreads) private(                  \
        a_data, wk_data, fft_data_data, Spower_data, Spower_iter_data, r2,     \
            b_a_data, b_wk_data, y_data, r3, r4, r5, y, r6, r7, b_y, Spower,   \
            wk, a, r8, fft_data, b_Spower, b_a, Spower_iter, b_wk, Tpower,     \
            emlrtJBEnviron, g_st, i1, c_window_start, b_loop_ub, d, i2,        \
            c_loop_ub, d_loop_ub, b_ii, e_loop_ub, f_loop_ub)                  \
    firstprivate(d_st, e_st, f_st, emlrtHadParallelError)
  {
    if (setjmp(emlrtJBEnviron) == 0) {
      d_st.prev = sp;
      d_st.tls = emlrtAllocTLS((emlrtCTX)sp, omp_get_thread_num());
      d_st.site = NULL;
      emlrtSetJmpBuf(&d_st, &emlrtJBEnviron);
      e_st.prev = &d_st;
      e_st.tls = d_st.tls;
      f_st.prev = &e_st;
      f_st.tls = e_st.tls;
      g_st.prev = &f_st;
      g_st.tls = f_st.tls;
      emxInit_int32_T(&d_st, &r5, 2, &rc_emlrtRTEI);
      emxInit_real32_T(&d_st, &y, 1, &qc_emlrtRTEI);
      emxInit_real32_T(&d_st, &r6, 2, &jc_emlrtRTEI);
      emxInit_real_T(&d_st, &r7, 2, &pc_emlrtRTEI);
      emxInit_real32_T(&d_st, &b_y, 2, &ec_emlrtRTEI);
      emxInit_real32_T(&d_st, &Spower, 2, &ac_emlrtRTEI);
      emxInit_real32_T(&d_st, &wk, 2, &ic_emlrtRTEI);
      emxInit_boolean_T(&d_st, &a, 1, &qb_emlrtRTEI);
      emxInit_int32_T(&d_st, &r8, 2, &oc_emlrtRTEI);
      r2 = r8->data;
      emxInit_creal32_T(&d_st, &fft_data, &nc_emlrtRTEI);
      emxInit_real32_T(&d_st, &b_Spower, 2, &mc_emlrtRTEI);
      emxInit_real32_T(&d_st, &b_a, 1, &cc_emlrtRTEI);
      emxInit_real32_T(&d_st, &Spower_iter, 1, &lc_emlrtRTEI);
      emxInit_real32_T(&d_st, &b_wk, 2, &kc_emlrtRTEI);
    } else {
      emlrtHadParallelError = true;
    }
#pragma omp for nowait
    for (n = 0; n < input_sizes_idx_0; n++) {
      if (emlrtHadParallelError) {
        continue;
      }
      if (setjmp(emlrtJBEnviron) == 0) {
        /* Grab the data for the given window */
        if (winsize_samples - 1 < 0) {
          r8->size[0] = 1;
          r8->size[1] = 0;
        } else {
          i1 = r8->size[0] * r8->size[1];
          r8->size[0] = 1;
          r8->size[1] = winsize_samples;
          emxEnsureCapacity_int32_T(&d_st, r8, i1, &mb_emlrtRTEI);
          r2 = r8->data;
          for (i1 = 0; i1 < winsize_samples; i1++) {
            r2[i1] = i1;
          }
        }
        e_st.site = &e_emlrtRSI;
        indexShapeCheck(&e_st, data->size[0], r8->size);
        if (n + 1 > window_start->size[1]) {
          emlrtDynamicBoundsCheckR2012b(n + 1, 1, window_start->size[1],
                                        &d_emlrtBCI, &d_st);
        }
        c_window_start = window_start_data[n];
        i1 = b_a->size[0];
        b_a->size[0] = r8->size[1];
        emxEnsureCapacity_real32_T(&d_st, b_a, i1, &ob_emlrtRTEI);
        a_data = b_a->data;
        b_loop_ub = r8->size[1];
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          d = c_window_start + (real_T)r2[i1];
          if (d != (int32_T)muDoubleScalarFloor(d)) {
            emlrtIntegerCheckR2012b(d, &d_emlrtDCI, &d_st);
          }
          if (((int32_T)d < 1) || ((int32_T)d > data->size[0])) {
            emlrtDynamicBoundsCheckR2012b((int32_T)d, 1, data->size[0],
                                          &e_emlrtBCI, &d_st);
          }
          a_data[i1] = data_data[(int32_T)d - 1];
        }
        /* Skip empty segments */
        i1 = a->size[0];
        a->size[0] = b_a->size[0];
        emxEnsureCapacity_boolean_T(&d_st, a, i1, &qb_emlrtRTEI);
        b_a_data = a->data;
        b_loop_ub = b_a->size[0];
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          b_a_data[i1] = (a_data[i1] == 0.0F);
        }
        e_st.site = &f_emlrtRSI;
        if (!all(&e_st, a)) {
          i1 = a->size[0];
          a->size[0] = b_a->size[0];
          emxEnsureCapacity_boolean_T(&d_st, a, i1, &rb_emlrtRTEI);
          b_a_data = a->data;
          b_loop_ub = b_a->size[0];
          for (i1 = 0; i1 < b_loop_ub; i1++) {
            b_a_data[i1] = muSingleScalarIsNaN(a_data[i1]);
          }
          e_st.site = &g_emlrtRSI;
          if (any(&e_st, a)) {
            if (n + 1 > mt_spectrogram->size[1]) {
              emlrtDynamicBoundsCheckR2012b(n + 1, 1, mt_spectrogram->size[1],
                                            &b_emlrtBCI, &d_st);
            }
            b_loop_ub = mt_spectrogram->size[0];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              mt_spectrogram_data[i1 + mt_spectrogram->size[0] * n] = rtNaNF;
            }
          } else {
            /* Option to detrend_opt data to remove low frequency DC component
             */
            if (detrend_opt == 1.0) {
              i1 = Spower_iter->size[0];
              Spower_iter->size[0] = b_a->size[0];
              emxEnsureCapacity_real32_T(&d_st, Spower_iter, i1, &tb_emlrtRTEI);
              Spower_iter_data = Spower_iter->data;
              b_loop_ub = b_a->size[0] - 1;
              for (i1 = 0; i1 <= b_loop_ub; i1++) {
                Spower_iter_data[i1] = a_data[i1];
              }
              e_st.site = &h_emlrtRSI;
              detrend(&e_st, Spower_iter, b_a);
            } else if (detrend_opt == 2.0) {
              e_st.site = &i_emlrtRSI;
              b_detrend(&e_st, b_a);
            }
            /* Multiply the data by the tapers (STEP 2) */
            e_st.site = &j_emlrtRSI;
            repmat(&e_st, b_a, num_tapers, b_wk);
            wk_data = b_wk->data;
            if ((b_wk->size[0] != DPSS_tapers->size[0]) &&
                ((b_wk->size[0] != 1) && (DPSS_tapers->size[0] != 1))) {
              emlrtDimSizeImpxCheckR2021b(b_wk->size[0], DPSS_tapers->size[0],
                                          &i_emlrtECI, &d_st);
            }
            if ((b_wk->size[1] != DPSS_tapers->size[1]) &&
                ((b_wk->size[1] != 1) && (DPSS_tapers->size[1] != 1))) {
              emlrtDimSizeImpxCheckR2021b(b_wk->size[1], DPSS_tapers->size[1],
                                          &h_emlrtECI, &d_st);
            }
            if ((b_wk->size[0] == DPSS_tapers->size[0]) &&
                (b_wk->size[1] == DPSS_tapers->size[1])) {
              b_loop_ub = b_wk->size[0] * b_wk->size[1];
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                wk_data[i1] *= (real32_T)DPSS_tapers_data[i1];
              }
            } else {
              e_st.site = &j_emlrtRSI;
              binary_expand_op_2(&e_st, b_wk, DPSS_tapers);
            }
            /* Compute the FFT (STEP 3) */
            e_st.site = &k_emlrtRSI;
            fft(&e_st, b_wk, nfft, fft_data);
            fft_data_data = fft_data->data;
            /* Compute the weighted mean spectral power across tapers (STEP 4)
             */
            e_st.site = &l_emlrtRSI;
            f_st.site = &cb_emlrtRSI;
            i1 = b_Spower->size[0] * b_Spower->size[1];
            b_Spower->size[0] = fft_data->size[0];
            b_Spower->size[1] = fft_data->size[1];
            emxEnsureCapacity_real32_T(&f_st, b_Spower, i1, &wb_emlrtRTEI);
            Spower_data = b_Spower->data;
            b_loop_ub = fft_data->size[0] * fft_data->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              Tpower = fft_data_data[i1].im;
              Spower_data[i1] = Tpower * Tpower;
            }
            e_st.site = &l_emlrtRSI;
            f_st.site = &cb_emlrtRSI;
            i1 = wk->size[0] * wk->size[1];
            wk->size[0] = fft_data->size[0];
            wk->size[1] = fft_data->size[1];
            emxEnsureCapacity_real32_T(&f_st, wk, i1, &xb_emlrtRTEI);
            b_wk_data = wk->data;
            b_loop_ub = fft_data->size[0] * fft_data->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              Tpower = fft_data_data[i1].re;
              b_wk_data[i1] = Tpower * Tpower;
            }
            if ((b_Spower->size[0] != wk->size[0]) ||
                (b_Spower->size[1] != wk->size[1])) {
              emlrtSizeEqCheckNDErrorR2021b(&b_Spower->size[0], &wk->size[0],
                                            &g_emlrtECI, &d_st);
            }
            b_loop_ub = b_Spower->size[0] * b_Spower->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              Spower_data[i1] += b_wk_data[i1];
            }
            if (weighting == 2.0) {
              /*  adaptive weights - for colored noise spectrum (Percival &
               * Walden */
              /*  p368-p370) */
              e_st.site = &m_emlrtRSI;
              f_st.site = &nd_emlrtRSI;
              dynamic_size_checks(&f_st, b_a, b_a, b_a->size[0], b_a->size[0]);
              Tpower = mtimes(b_a, b_a) / (real32_T)b_a->size[0];
              i1 = Spower->size[0] * Spower->size[1];
              Spower->size[0] = b_Spower->size[0];
              Spower->size[1] = 2;
              emxEnsureCapacity_real32_T(&d_st, Spower, i1, &ac_emlrtRTEI);
              wk_data = Spower->data;
              b_loop_ub = b_Spower->size[0];
              for (i1 = 0; i1 < 2; i1++) {
                for (i2 = 0; i2 < b_loop_ub; i2++) {
                  if (i1 + 1 > b_Spower->size[1]) {
                    emlrtDynamicBoundsCheckR2012b(i1 + 1, 1, b_Spower->size[1],
                                                  &i_emlrtBCI, &d_st);
                  }
                  wk_data[i2 + Spower->size[0] * i1] =
                      Spower_data[i2 + b_Spower->size[0] * i1];
                }
              }
              e_st.site = &n_emlrtRSI;
              mean(&e_st, Spower, Spower_iter);
              Spower_iter_data = Spower_iter->data;
              b_loop_ub = DPSS_eigen->size[0];
              i1 = b_a->size[0];
              b_a->size[0] = DPSS_eigen->size[0];
              emxEnsureCapacity_real32_T(&d_st, b_a, i1, &cc_emlrtRTEI);
              a_data = b_a->data;
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                a_data[i1] = (real32_T)(1.0 - DPSS_eigen_data[i1]) * Tpower;
              }
              b_loop_ub = DPSS_eigen->size[0];
              c_loop_ub = b_a->size[0];
              d = (int32_T)muDoubleScalarFloor(nfft);
              d_loop_ub = (int32_T)nfft;
              for (b_ii = 0; b_ii < 3; b_ii++) {
                /*  run 3 iterations */
                /*  calculate the MSE weights */
                i1 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = Spower_iter->size[0];
                b_y->size[1] = DPSS_eigen->size[0];
                emxEnsureCapacity_real32_T(&d_st, b_y, i1, &ec_emlrtRTEI);
                y_data = b_y->data;
                for (i1 = 0; i1 < b_loop_ub; i1++) {
                  e_loop_ub = Spower_iter->size[0];
                  for (i2 = 0; i2 < e_loop_ub; i2++) {
                    y_data[i2 + b_y->size[0] * i1] =
                        Spower_iter_data[i2] * (real32_T)DPSS_eigen_data[i1];
                  }
                }
                if (nfft != d) {
                  emlrtIntegerCheckR2012b(nfft, &b_emlrtDCI, &d_st);
                }
                i1 = wk->size[0] * wk->size[1];
                wk->size[0] = (int32_T)nfft;
                wk->size[1] = b_a->size[0];
                emxEnsureCapacity_real32_T(&d_st, wk, i1, &fc_emlrtRTEI);
                b_wk_data = wk->data;
                for (i1 = 0; i1 < c_loop_ub; i1++) {
                  for (i2 = 0; i2 < d_loop_ub; i2++) {
                    b_wk_data[i2 + wk->size[0] * i1] = a_data[i1];
                  }
                }
                if ((b_y->size[0] != wk->size[0]) &&
                    ((b_y->size[0] != 1) && (wk->size[0] != 1))) {
                  emlrtDimSizeImpxCheckR2021b(b_y->size[0], wk->size[0],
                                              &f_emlrtECI, &d_st);
                }
                i1 = b_wk->size[0] * b_wk->size[1];
                b_wk->size[0] = Spower_iter->size[0];
                b_wk->size[1] = num_tapers;
                emxEnsureCapacity_real32_T(&d_st, b_wk, i1, &gc_emlrtRTEI);
                wk_data = b_wk->data;
                for (i1 = 0; i1 < num_tapers; i1++) {
                  e_loop_ub = Spower_iter->size[0];
                  for (i2 = 0; i2 < e_loop_ub; i2++) {
                    wk_data[i2 + b_wk->size[0] * i1] = Spower_iter_data[i2];
                  }
                }
                e_st.site = &o_emlrtRSI;
                if (b_y->size[0] == wk->size[0]) {
                  e_loop_ub = b_y->size[0] * b_y->size[1];
                  for (i1 = 0; i1 < e_loop_ub; i1++) {
                    y_data[i1] += b_wk_data[i1];
                  }
                } else {
                  f_st.site = &o_emlrtRSI;
                  plus(&f_st, b_y, wk);
                  y_data = b_y->data;
                }
                f_st.site = &td_emlrtRSI;
                g_st.site = &ud_emlrtRSI;
                assertCompatibleDims(&g_st, b_wk, b_y);
                if ((b_wk->size[0] == b_y->size[0]) &&
                    (b_wk->size[1] == b_y->size[1])) {
                  e_loop_ub = b_wk->size[0] * b_wk->size[1];
                  for (i1 = 0; i1 < e_loop_ub; i1++) {
                    wk_data[i1] /= y_data[i1];
                  }
                } else {
                  g_st.site = &bf_emlrtRSI;
                  b_rdivide(&g_st, b_wk, b_y);
                  wk_data = b_wk->data;
                }
                /*  calculate new spectral estimate */
                e_st.site = &p_emlrtRSI;
                f_st.site = &cb_emlrtRSI;
                e_loop_ub = b_wk->size[0] * b_wk->size[1];
                for (i1 = 0; i1 < e_loop_ub; i1++) {
                  Tpower = wk_data[i1];
                  wk_data[i1] = Tpower * Tpower;
                }
                if (nfft != d) {
                  emlrtIntegerCheckR2012b(nfft, &emlrtDCI, &d_st);
                }
                i1 = r7->size[0] * r7->size[1];
                r7->size[0] = (int32_T)nfft;
                r7->size[1] = DPSS_eigen->size[0];
                emxEnsureCapacity_real_T(&d_st, r7, i1, &hc_emlrtRTEI);
                r3 = r7->data;
                for (i1 = 0; i1 < b_loop_ub; i1++) {
                  for (i2 = 0; i2 < d_loop_ub; i2++) {
                    r3[i2 + r7->size[0] * i1] = DPSS_eigen_data[i1];
                  }
                }
                if ((b_wk->size[0] != r7->size[0]) &&
                    ((b_wk->size[0] != 1) && (r7->size[0] != 1))) {
                  emlrtDimSizeImpxCheckR2021b(b_wk->size[0], r7->size[0],
                                              &e_emlrtECI, &d_st);
                }
                if ((b_wk->size[1] != r7->size[1]) &&
                    ((b_wk->size[1] != 1) && (r7->size[1] != 1))) {
                  emlrtDimSizeImpxCheckR2021b(b_wk->size[1], r7->size[1],
                                              &d_emlrtECI, &d_st);
                }
                if ((b_wk->size[0] == r7->size[0]) &&
                    (b_wk->size[1] == r7->size[1])) {
                  e_loop_ub = b_wk->size[0] * b_wk->size[1];
                  for (i1 = 0; i1 < e_loop_ub; i1++) {
                    wk_data[i1] *= (real32_T)r3[i1];
                  }
                } else {
                  e_st.site = &p_emlrtRSI;
                  binary_expand_op_1(&e_st, b_wk, r7);
                  wk_data = b_wk->data;
                }
                if ((b_wk->size[1] != b_Spower->size[1]) &&
                    ((b_wk->size[1] != 1) && (b_Spower->size[1] != 1))) {
                  emlrtDimSizeImpxCheckR2021b(b_wk->size[1], b_Spower->size[1],
                                              &c_emlrtECI, &d_st);
                }
                if ((b_wk->size[0] != b_Spower->size[0]) &&
                    ((b_wk->size[0] != 1) && (b_Spower->size[0] != 1))) {
                  emlrtDimSizeImpxCheckR2021b(b_wk->size[0], b_Spower->size[0],
                                              &b_emlrtECI, &d_st);
                }
                e_st.site = &q_emlrtRSI;
                if ((b_wk->size[1] == b_Spower->size[1]) &&
                    (b_wk->size[0] == b_Spower->size[0])) {
                  i1 = wk->size[0] * wk->size[1];
                  wk->size[0] = b_wk->size[1];
                  wk->size[1] = b_wk->size[0];
                  emxEnsureCapacity_real32_T(&e_st, wk, i1, &ic_emlrtRTEI);
                  b_wk_data = wk->data;
                  e_loop_ub = b_wk->size[0];
                  for (i1 = 0; i1 < e_loop_ub; i1++) {
                    f_loop_ub = b_wk->size[1];
                    for (i2 = 0; i2 < f_loop_ub; i2++) {
                      b_wk_data[i2 + wk->size[0] * i1] =
                          wk_data[i1 + b_wk->size[0] * i2] *
                          Spower_data[i1 + b_Spower->size[0] * i2];
                    }
                  }
                  f_st.site = &q_emlrtRSI;
                  sum(&f_st, wk, r6);
                  wk_data = r6->data;
                } else {
                  f_st.site = &q_emlrtRSI;
                  binary_expand_op(&f_st, r6, q_emlrtRSI, b_wk, b_Spower);
                  wk_data = r6->data;
                }
                i1 = Spower_iter->size[0];
                Spower_iter->size[0] = r6->size[1];
                emxEnsureCapacity_real32_T(&e_st, Spower_iter, i1,
                                           &jc_emlrtRTEI);
                Spower_iter_data = Spower_iter->data;
                e_loop_ub = r6->size[1];
                for (i1 = 0; i1 < e_loop_ub; i1++) {
                  Spower_iter_data[i1] = wk_data[i1];
                }
                f_st.site = &q_emlrtRSI;
                b_sum(&f_st, b_wk, y);
                y_data = y->data;
                f_st.site = &td_emlrtRSI;
                g_st.site = &ud_emlrtRSI;
                b_assertCompatibleDims(&g_st, Spower_iter, y);
                e_loop_ub = Spower_iter->size[0];
                if (Spower_iter->size[0] == y->size[0]) {
                  for (i1 = 0; i1 < e_loop_ub; i1++) {
                    Spower_iter_data[i1] /= y_data[i1];
                  }
                } else {
                  g_st.site = &bf_emlrtRSI;
                  rdivide(&g_st, Spower_iter, y);
                  Spower_iter_data = Spower_iter->data;
                }
                if (*emlrtBreakCheckR2012bFlagVar != 0) {
                  emlrtBreakCheckR2012b(&d_st);
                }
              }
            } else {
              /*  eigenvalue or uniform weights */
              e_st.site = &r_emlrtRSI;
              f_st.site = &nd_emlrtRSI;
              b_dynamic_size_checks(&f_st, b_Spower, wt, b_Spower->size[1],
                                    wt->size[0]);
              i1 = Spower_iter->size[0];
              Spower_iter->size[0] = b_Spower->size[0];
              emxEnsureCapacity_real32_T(&e_st, Spower_iter, i1, &yb_emlrtRTEI);
              Spower_iter_data = Spower_iter->data;
              b_loop_ub = b_Spower->size[0];
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                Spower_iter_data[i1] = 0.0F;
                c_loop_ub = b_Spower->size[1];
                for (i2 = 0; i2 < c_loop_ub; i2++) {
                  Spower_iter_data[i1] +=
                      Spower_data[i1 + b_Spower->size[0] * i2] *
                      (real32_T)stimes_data[i2];
                }
              }
            }
            /* Append the spectrum to the spectrogram */
            if (n + 1 > mt_spectrogram->size[1]) {
              emlrtDynamicBoundsCheckR2012b(n + 1, 1, mt_spectrogram->size[1],
                                            &emlrtBCI, &d_st);
            }
            c_loop_ub = freq_inds->size[1];
            b_loop_ub = 0;
            for (d_loop_ub = 0; d_loop_ub < c_loop_ub; d_loop_ub++) {
              if (freq_inds_data[d_loop_ub]) {
                b_loop_ub++;
              }
            }
            i1 = r5->size[0] * r5->size[1];
            r5->size[0] = 1;
            r5->size[1] = b_loop_ub;
            emxEnsureCapacity_int32_T(&d_st, r5, i1, &db_emlrtRTEI);
            r4 = r5->data;
            b_loop_ub = 0;
            for (d_loop_ub = 0; d_loop_ub < c_loop_ub; d_loop_ub++) {
              if (freq_inds_data[d_loop_ub]) {
                r4[b_loop_ub] = d_loop_ub;
                b_loop_ub++;
              }
            }
            b_loop_ub = r5->size[1];
            i1 = b_a->size[0];
            b_a->size[0] = r5->size[1];
            emxEnsureCapacity_real32_T(&d_st, b_a, i1, &dc_emlrtRTEI);
            a_data = b_a->data;
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              if ((r4[i1] < 0) || (r4[i1] > Spower_iter->size[0] - 1)) {
                emlrtDynamicBoundsCheckR2012b(
                    r4[i1], 0, Spower_iter->size[0] - 1, &j_emlrtBCI, &d_st);
              }
              a_data[i1] = Spower_iter_data[r4[i1]];
            }
            emlrtSubAssignSizeCheckR2012b(&mt_spectrogram->size[0], 1,
                                          &r5->size[1], 1, &emlrtECI, &d_st);
            b_loop_ub = mt_spectrogram->size[0];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              mt_spectrogram_data[i1 + mt_spectrogram->size[0] * n] =
                  a_data[i1];
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
      emlrtHeapReferenceStackLeaveScope(&d_st, 14);
      emxFree_real32_T(&d_st, &b_wk);
      emxFree_real32_T(&d_st, &Spower_iter);
      emxFree_real32_T(&d_st, &b_a);
      emxFree_real32_T(&d_st, &b_Spower);
      emxFree_creal32_T(&d_st, &fft_data);
      emxFree_int32_T(&d_st, &r8);
      emxFree_boolean_T(&d_st, &a);
      emxFree_real32_T(&d_st, &wk);
      emxFree_real32_T(&d_st, &Spower);
      emxFree_real32_T(&d_st, &b_y);
      emxFree_real_T(&d_st, &r7);
      emxFree_real32_T(&d_st, &r6);
      emxFree_real32_T(&d_st, &y);
      emxFree_int32_T(&d_st, &r5);
    }
  }
  emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
  emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  emxFree_real_T(sp, &wt);
  /* Compute one-sided PSD spectrum  */
  st.site = &s_emlrtRSI;
  i = freq_inds->size[0] * freq_inds->size[1];
  freq_inds->size[0] = 1;
  sizes_idx_0 = sfreqs->size[1];
  freq_inds->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(&st, freq_inds, i, &jb_emlrtRTEI);
  freq_inds_data = freq_inds->data;
  for (i = 0; i < sizes_idx_0; i++) {
    freq_inds_data[i] = (sfreqs_data[i] == 0.0);
  }
  emxInit_int32_T(&st, &ii, 2, &uc_emlrtRTEI);
  b_st.site = &oe_emlrtRSI;
  eml_find(&b_st, freq_inds, ii);
  ii_data = ii->data;
  i = window_start->size[0] * window_start->size[1];
  window_start->size[0] = 1;
  b_input_sizes_idx_0 = ii->size[1];
  window_start->size[1] = ii->size[1];
  emxEnsureCapacity_real_T(&st, window_start, i, &kb_emlrtRTEI);
  window_start_data = window_start->data;
  for (i = 0; i < b_input_sizes_idx_0; i++) {
    window_start_data[i] = ii_data[i];
  }
  st.site = &t_emlrtRSI;
  b = Fs / 2.0;
  i = freq_inds->size[0] * freq_inds->size[1];
  freq_inds->size[0] = 1;
  freq_inds->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(&st, freq_inds, i, &lb_emlrtRTEI);
  freq_inds_data = freq_inds->data;
  for (i = 0; i < sizes_idx_0; i++) {
    freq_inds_data[i] = (sfreqs_data[i] == b);
  }
  b_st.site = &oe_emlrtRSI;
  eml_find(&b_st, freq_inds, ii);
  ii_data = ii->data;
  emxFree_boolean_T(&st, &freq_inds);
  emxInit_real_T(sp, &c_y, 2, &vc_emlrtRTEI);
  if (sfreqs->size[1] < 1) {
    c_y->size[0] = 1;
    c_y->size[1] = 0;
  } else {
    i = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = sfreqs->size[1];
    emxEnsureCapacity_real_T(sp, c_y, i, &nb_emlrtRTEI);
    stimes_data = c_y->data;
    sizes_idx_0 = sfreqs->size[1] - 1;
    for (i = 0; i <= sizes_idx_0; i++) {
      stimes_data[i] = (real_T)i + 1.0;
    }
  }
  st.site = &u_emlrtRSI;
  b_st.site = &re_emlrtRSI;
  emxInit_real_T(&b_st, &b_window_start, 2, &pb_emlrtRTEI);
  i = b_window_start->size[0] * b_window_start->size[1];
  b_window_start->size[0] = 1;
  b_window_start->size[1] = window_start->size[1] + ii->size[1];
  emxEnsureCapacity_real_T(&b_st, b_window_start, i, &pb_emlrtRTEI);
  stimes_data = b_window_start->data;
  for (i = 0; i < b_input_sizes_idx_0; i++) {
    stimes_data[i] = window_start_data[i];
  }
  b_i = ii->size[1];
  for (i = 0; i < b_i; i++) {
    stimes_data[i + window_start->size[1]] = ii_data[i];
  }
  emxInit_real_T(&b_st, &b_select, 2, &sc_emlrtRTEI);
  emxInit_int32_T(&b_st, &ia, 1, &db_emlrtRTEI);
  c_st.site = &se_emlrtRSI;
  do_vectors(&c_st, c_y, b_window_start, b_select, ia);
  stimes_data = b_select->data;
  emxFree_real_T(&b_st, &b_window_start);
  emxFree_int32_T(&b_st, &ia);
  emxFree_real_T(&b_st, &c_y);
  emxInit_real32_T(sp, &b_b, 2, &sb_emlrtRTEI);
  sizes_idx_0 = b_select->size[1];
  i = b_b->size[0] * b_b->size[1];
  b_b->size[0] = b_select->size[1];
  b_b->size[1] = loop_ub;
  emxEnsureCapacity_real32_T(sp, b_b, i, &sb_emlrtRTEI);
  b_data = b_b->data;
  for (i = 0; i < loop_ub; i++) {
    for (winsize_samples = 0; winsize_samples < sizes_idx_0;
         winsize_samples++) {
      if (stimes_data[winsize_samples] !=
          (int32_T)muDoubleScalarFloor(stimes_data[winsize_samples])) {
        emlrtIntegerCheckR2012b(stimes_data[winsize_samples], &e_emlrtDCI,
                                (emlrtConstCTX)sp);
      }
      input_sizes_idx_0 = (int32_T)stimes_data[winsize_samples];
      if ((input_sizes_idx_0 < 1) || (input_sizes_idx_0 > nz_tmp)) {
        emlrtDynamicBoundsCheckR2012b(input_sizes_idx_0, 1, nz_tmp, &f_emlrtBCI,
                                      (emlrtConstCTX)sp);
      }
      b_data[winsize_samples + b_b->size[0] * i] = mt_spectrogram_data
          [(input_sizes_idx_0 + mt_spectrogram->size[0] * i) - 1];
    }
  }
  emxFree_real_T(sp, &b_select);
  input_sizes_idx_0 = b_b->size[0] * b_b->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    b_data[i] *= 2.0F;
  }
  st.site = &v_emlrtRSI;
  emxInit_real32_T(&st, &varargin_1, 2, &ub_emlrtRTEI);
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = b_input_sizes_idx_0;
  varargin_1->size[1] = loop_ub;
  emxEnsureCapacity_real32_T(&st, varargin_1, i, &ub_emlrtRTEI);
  varargin_1_data = varargin_1->data;
  for (i = 0; i < loop_ub; i++) {
    for (winsize_samples = 0; winsize_samples < b_input_sizes_idx_0;
         winsize_samples++) {
      input_sizes_idx_0 = (int32_T)window_start_data[winsize_samples];
      if ((input_sizes_idx_0 < 1) || (input_sizes_idx_0 > nz_tmp)) {
        emlrtDynamicBoundsCheckR2012b(input_sizes_idx_0, 1, nz_tmp, &g_emlrtBCI,
                                      &st);
      }
      varargin_1_data[winsize_samples + varargin_1->size[0] * i] =
          mt_spectrogram_data[(input_sizes_idx_0 +
                               mt_spectrogram->size[0] * i) -
                              1];
    }
  }
  emxInit_real32_T(&st, &varargin_3, 2, &vb_emlrtRTEI);
  i = varargin_3->size[0] * varargin_3->size[1];
  varargin_3->size[0] = ii->size[1];
  varargin_3->size[1] = loop_ub;
  emxEnsureCapacity_real32_T(&st, varargin_3, i, &vb_emlrtRTEI);
  varargin_3_data = varargin_3->data;
  for (i = 0; i < loop_ub; i++) {
    for (winsize_samples = 0; winsize_samples < b_i; winsize_samples++) {
      if ((ii_data[winsize_samples] < 1) ||
          (ii_data[winsize_samples] > nz_tmp)) {
        emlrtDynamicBoundsCheckR2012b(ii_data[winsize_samples], 1, nz_tmp,
                                      &h_emlrtBCI, &st);
      }
      varargin_3_data[winsize_samples + varargin_3->size[0] * i] =
          mt_spectrogram_data[(ii_data[winsize_samples] +
                               mt_spectrogram->size[0] * i) -
                              1];
    }
  }
  b_st.site = &xe_emlrtRSI;
  c_b = ((window_start->size[1] != 0) && (mt_spectrogram->size[1] != 0));
  if (c_b) {
    b_i = loop_ub;
  } else if ((b_b->size[0] != 0) && (b_b->size[1] != 0)) {
    b_i = b_b->size[1];
  } else if ((ii->size[1] != 0) && (mt_spectrogram->size[1] != 0)) {
    b_i = loop_ub;
  } else {
    b_i = loop_ub;
    if (b_b->size[1] > loop_ub) {
      b_i = b_b->size[1];
    }
    if (mt_spectrogram->size[1] > b_i) {
      b_i = loop_ub;
    }
  }
  c_st.site = &ye_emlrtRSI;
  if ((mt_spectrogram->size[1] != b_i) &&
      ((window_start->size[1] != 0) && (mt_spectrogram->size[1] != 0))) {
    emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  if ((b_b->size[1] != b_i) && ((b_b->size[0] != 0) && (b_b->size[1] != 0))) {
    emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  if ((mt_spectrogram->size[1] != b_i) &&
      ((ii->size[1] != 0) && (mt_spectrogram->size[1] != 0))) {
    emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  empty_non_axis_sizes = (b_i == 0);
  if (empty_non_axis_sizes || c_b) {
    input_sizes_idx_0 = window_start->size[1];
  } else {
    input_sizes_idx_0 = 0;
  }
  emxFree_real_T(&b_st, &window_start);
  if (empty_non_axis_sizes || ((b_b->size[0] != 0) && (b_b->size[1] != 0))) {
    b_input_sizes_idx_0 = b_b->size[0];
  } else {
    b_input_sizes_idx_0 = 0;
  }
  if (empty_non_axis_sizes ||
      ((ii->size[1] != 0) && (mt_spectrogram->size[1] != 0))) {
    sizes_idx_0 = ii->size[1];
  } else {
    sizes_idx_0 = 0;
  }
  emxFree_int32_T(&b_st, &ii);
  num_tapers_tmp = input_sizes_idx_0;
  i = mt_spectrogram->size[0] * mt_spectrogram->size[1];
  mt_spectrogram->size[0] =
      (num_tapers_tmp + b_input_sizes_idx_0) + sizes_idx_0;
  mt_spectrogram->size[1] = b_i;
  emxEnsureCapacity_real32_T(sp, mt_spectrogram, i, &bc_emlrtRTEI);
  mt_spectrogram_data = mt_spectrogram->data;
  for (i = 0; i < b_i; i++) {
    for (winsize_samples = 0; winsize_samples < num_tapers_tmp;
         winsize_samples++) {
      mt_spectrogram_data[winsize_samples + mt_spectrogram->size[0] * i] =
          varargin_1_data[winsize_samples + num_tapers_tmp * i] / (real32_T)Fs;
    }
    for (winsize_samples = 0; winsize_samples < b_input_sizes_idx_0;
         winsize_samples++) {
      mt_spectrogram_data[(winsize_samples + num_tapers_tmp) +
                          mt_spectrogram->size[0] * i] =
          b_data[winsize_samples + b_input_sizes_idx_0 * i] / (real32_T)Fs;
    }
    for (winsize_samples = 0; winsize_samples < sizes_idx_0;
         winsize_samples++) {
      mt_spectrogram_data[((winsize_samples + num_tapers_tmp) +
                           b_input_sizes_idx_0) +
                          mt_spectrogram->size[0] * i] =
          varargin_3_data[winsize_samples + sizes_idx_0 * i] / (real32_T)Fs;
    }
  }
  emxFree_real32_T(sp, &varargin_3);
  emxFree_real32_T(sp, &varargin_1);
  emxFree_real32_T(sp, &b_b);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (multitaper_spectrogram_coder.c) */
