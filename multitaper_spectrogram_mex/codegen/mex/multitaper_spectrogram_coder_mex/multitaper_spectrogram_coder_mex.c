/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * multitaper_spectrogram_coder_mex.c
 *
 * Code generation for function 'multitaper_spectrogram_coder_mex'
 *
 */

/* Include files */
#include "multitaper_spectrogram_coder_mex.h"
#include "blockedSummation.h"
#include "colon.h"
#include "combineVectorElements.h"
#include "detrend.h"
#include "eml_int_forloop_overflow_check.h"
#include "fft.h"
#include "indexShapeCheck.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex_emxutil.h"
#include "mwmathutil.h"
#include "nextpow2.h"
#include "power.h"
#include "repmat.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 69,    /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 74,  /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 78,  /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 95,  /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 110, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 113, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 117, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 125, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 127, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 131, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 134, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 138, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 139, /* lineNo */
  "multitaper_spectrogram_coder_mex",  /* fcnName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 103, /* lineNo */
  "colon",                             /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/ops/colon.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 15,  /* lineNo */
  "sum",                               /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 96, /* lineNo */
  "sumprod",                           /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pathName */
};

static emlrtRSInfo eb_emlrtRSI = { 13, /* lineNo */
  "all",                               /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/ops/all.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 143,/* lineNo */
  "allOrAny",                          /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/allOrAny.m"/* pathName */
};

static emlrtRSInfo gb_emlrtRSI = { 13, /* lineNo */
  "any",                               /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/ops/any.m"/* pathName */
};

static emlrtRSInfo we_emlrtRSI = { 20, /* lineNo */
  "sum",                               /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

static emlrtRSInfo xe_emlrtRSI = { 62, /* lineNo */
  "combineVectorElements",             /* fcnName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  135,                                 /* lineNo */
  26,                                  /* colNo */
  "fft_data",                          /* aName */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 110,   /* lineNo */
  25,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  110,                                 /* lineNo */
  25,                                  /* colNo */
  "data",                              /* aName */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  110,                                 /* lineNo */
  25,                                  /* colNo */
  "window_start",                      /* aName */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  87,                                  /* lineNo */
  10,                                  /* colNo */
  "sfreqs",                            /* aName */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 95,  /* lineNo */
  24,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m",/* pName */
  4                                    /* checkKind */
};

static emlrtECInfo emlrtECI = { 2,     /* nDims */
  86,                                  /* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtECInfo b_emlrtECI = { 2,   /* nDims */
  131,                                 /* lineNo */
  20,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtECInfo c_emlrtECI = { 2,   /* nDims */
  138,                                 /* lineNo */
  17,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  118,                                 /* lineNo */
  26,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  142,                                 /* lineNo */
  22,                                  /* colNo */
  "mt_spectrogram",                    /* aName */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo d_emlrtECI = { -1,  /* nDims */
  142,                                 /* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo r_emlrtRTEI = { 69,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo s_emlrtRTEI = { 78,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo t_emlrtRTEI = { 86,/* lineNo */
  14,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo u_emlrtRTEI = { 86,/* lineNo */
  47,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo v_emlrtRTEI = { 86,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo w_emlrtRTEI = { 1, /* lineNo */
  43,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo x_emlrtRTEI = { 91,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo y_emlrtRTEI = { 86,/* lineNo */
  13,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = { 95,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = { 47,/* lineNo */
  9,                                   /* colNo */
  "div",                               /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/eml/+coder/+internal/div.m"/* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = { 146,/* lineNo */
  1,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo db_emlrtRTEI = { 110,/* lineNo */
  44,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = { 110,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = { 113,/* lineNo */
  12,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = { 17,/* lineNo */
  9,                                   /* colNo */
  "isnan",                             /* fName */
  "/apps/source/matlab/2019b/toolbox/eml/lib/matlab/elmat/isnan.m"/* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = { 125,/* lineNo */
  32,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = { 127,/* lineNo */
  32,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo jb_emlrtRTEI = { 135,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo kb_emlrtRTEI = { 138,/* lineNo */
  17,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = { 138,/* lineNo */
  38,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo mb_emlrtRTEI = { 139,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo nb_emlrtRTEI = { 138,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = { 134,/* lineNo */
  5,                                   /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = { 135,/* lineNo */
  26,                                  /* colNo */
  "multitaper_spectrogram_coder_mex",  /* fName */
  "/data/preraugp/code/labcode/multitaper/multitaper_spectrogram_mex/multitaper_spectrogram_coder_mex.m"/* pName */
};

/* Function Definitions */
void multitaper_spectrogram_coder_mex(const emlrtStack *sp, const
  emxArray_real32_T *data, real_T Fs, real_T frequency_range[2], const
  emxArray_real_T *DPSS_tapers, real_T time_bandwidth, real_T num_tapers, real_T
  winsize_samples, real_T winstep_samples, real_T min_NFFT, real_T detrend_opt,
  emxArray_real32_T *mt_spectrogram, emxArray_real_T *stimes, emxArray_real_T
  *sfreqs)
{
  int32_T nz;
  real_T b;
  emxArray_real_T *window_start;
  int32_T i;
  real_T nfft;
  emxArray_boolean_T *x;
  emxArray_boolean_T *r;
  emxArray_boolean_T *freq_inds;
  int32_T trueCount;
  int32_T b_i;
  int32_T partialTrueCount;
  int32_T n;
  emxArray_real32_T *r1;
  emxArray_real32_T *fft_range;
  emxArray_real32_T *b_mt_spectrogram;
  emxArray_int32_T *r2;
  emxArray_real32_T *data_segment;
  emxArray_boolean_T *b_x;
  emxArray_real_T *r3;
  emxArray_real32_T *b_data_segment;
  emxArray_creal32_T *fft_data;
  emxArray_creal32_T *b_fft_range;
  emxArray_real32_T *magnitude;
  emxArray_real32_T *mt_spectrum;
  int32_T i1;
  int32_T loop_ub;
  real_T b_window_start;
  int32_T b_data;
  real_T d;
  int32_T end;
  boolean_T y;
  boolean_T exitg1;
  int32_T c_window_start[1];
  emxArray_real32_T c_data_segment;
  int32_T d_window_start[1];
  jmp_buf * volatile emlrtJBStack;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  jmp_buf b_emlrtJBEnviron;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  boolean_T emlrtHadParallelError = false;
  (void)time_bandwidth;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /* MULTITAPER_SPECTROGRAM  Compute the multitaper spectrogram for time series data */
  /*  */
  /*    Usage: */
  /*    Direct input: */
  /*    [spect,stimes,sfreqs] = multitaper_spectrogram_coder(data, Fs, frequency_range, taper_params, window_params, min_NFFT, detrend_opt) */
  /*  */
  /*    Input: */
  /*    data: 1 x <number of samples> vector - time series data-- required */
  /*    Fs: double - sampling frequency in Hz  -- required */
  /*    frequency_range: 1x2 vector - [<min frequency>, <max frequency>] (default: [0 nyquist]) */
  /*    taper_params: 1x2 vector - [<time-halfbandwidth product>, <number of tapers>] (default: [5 9]) */
  /*    window_params: 1x2 vector - [window size (seconds), step size (seconds)] (default: [5 1]) */
  /*    detrend_opt: string - detrend data window ('linear' (default), 'constant', 'off'); */
  /*    min_NFFT: double - minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0) */
  /*    plot_on: boolean to plot results (default: true) */
  /*    verbose: boolean to display spectrogram properties (default: true) */
  /*  */
  /*    Output: */
  /*    spect: TxF matrix of spectral power */
  /*    stimes: 1XT vector of times for the center of the spectral bins */
  /*    sfreqs: 1XF vector of frequency bins for the spectrogram */
  /*  */
  /*    Example: */
  /*       Fs=200; %Sampling Frequency */
  /*       frequency_range=[0 25]; %Limit frequencies from .5 to 25 Hz */
  /*       taper_params=[3 5]; %Time bandwidth and number of tapers */
  /*       window_params=[4 1]; %Window size is 4s with step size of 1s */
  /*  */
  /*       %Generate sample chirp data */
  /*       t=1/Fs:1/Fs:600; %Create 10 minutes of data */
  /*       f_start=1;f_end=20; % Set chirp range in Hz */
  /*       data=chirp(t,f_start,t(end),f_end,'logarithmic'); */
  /*  */
  /*       %Compute the multitaper spectrogram */
  /*       [spect,stimes,sfreqs] = multitaper_spectrogram(data,Fs,frequency_range, taper_params, window_params); */
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
  /*    Last modified 1/11/2019 */
  /*  ******************************************************************** */
  /*  PROCESS DATA AND PARAMETERS */
  /* Fix error in frequency range */
  if (frequency_range[1] > Fs / 2.0) {
    frequency_range[1] = Fs / 2.0;
  }

  /* Generate DPSS tapers (STEP 1) */
  /*  DPSS_tapers = coder.const(@dpss, time_bandwidth, num_tapers) * sqrt(Fs); */
  /* Total data length */
  if ((data->size[0] == 0) || (data->size[1] == 0)) {
    nz = 0;
  } else {
    nz = muIntScalarMax_sint32(data->size[0], data->size[1]);
  }

  /* Window start indices */
  st.site = &emlrtRSI;
  b = ((real_T)nz - winsize_samples) + 1.0;
  emxInit_real_T(&st, &window_start, 2, &r_emlrtRTEI, true);
  if (muDoubleScalarIsNaN(winstep_samples) || muDoubleScalarIsNaN(b)) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &r_emlrtRTEI);
    window_start->data[0] = rtNaN;
  } else if ((winstep_samples == 0.0) || ((1.0 < b) && (winstep_samples < 0.0)) ||
             ((b < 1.0) && (winstep_samples > 0.0))) {
    window_start->size[0] = 1;
    window_start->size[1] = 0;
  } else if (muDoubleScalarIsInf(b) && (muDoubleScalarIsInf(winstep_samples) ||
              (1.0 == b))) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &r_emlrtRTEI);
    window_start->data[0] = rtNaN;
  } else if (muDoubleScalarIsInf(winstep_samples)) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    window_start->size[1] = 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &r_emlrtRTEI);
    window_start->data[0] = 1.0;
  } else if (muDoubleScalarFloor(winstep_samples) == winstep_samples) {
    i = window_start->size[0] * window_start->size[1];
    window_start->size[0] = 1;
    nz = (int32_T)muDoubleScalarFloor((b - 1.0) / winstep_samples);
    window_start->size[1] = nz + 1;
    emxEnsureCapacity_real_T(&st, window_start, i, &r_emlrtRTEI);
    for (i = 0; i <= nz; i++) {
      window_start->data[i] = winstep_samples * (real_T)i + 1.0;
    }
  } else {
    b_st.site = &n_emlrtRSI;
    eml_float_colon(&b_st, winstep_samples, b, window_start);
  }

  /* Number of windows */
  /* Number of points in the FFT */
  st.site = &b_emlrtRSI;
  b_st.site = &b_emlrtRSI;
  b = nextpow2(winsize_samples);
  b_st.site = &w_emlrtRSI;
  st.site = &b_emlrtRSI;
  b_st.site = &b_emlrtRSI;
  nfft = nextpow2(min_NFFT);
  b_st.site = &w_emlrtRSI;
  nfft = muDoubleScalarMax(muDoubleScalarMax(muDoubleScalarPower(2.0, b),
    winsize_samples), muDoubleScalarPower(2.0, nfft));

  /* Create the frequency vector */
  b = Fs / nfft;
  st.site = &c_emlrtRSI;
  if (muDoubleScalarIsNaN(b) || muDoubleScalarIsNaN(Fs)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &s_emlrtRTEI);
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
    emxEnsureCapacity_real_T(&st, sfreqs, i, &s_emlrtRTEI);
    sfreqs->data[0] = rtNaN;
  } else if (muDoubleScalarIsInf(b)) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    sfreqs->size[1] = 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &s_emlrtRTEI);
    sfreqs->data[0] = 0.0;
  } else if (muDoubleScalarFloor(b) == b) {
    i = sfreqs->size[0] * sfreqs->size[1];
    sfreqs->size[0] = 1;
    nz = (int32_T)muDoubleScalarFloor(Fs / b);
    sfreqs->size[1] = nz + 1;
    emxEnsureCapacity_real_T(&st, sfreqs, i, &s_emlrtRTEI);
    for (i = 0; i <= nz; i++) {
      sfreqs->data[i] = b * (real_T)i;
    }
  } else {
    b_st.site = &n_emlrtRSI;
    b_eml_float_colon(&b_st, b, Fs, sfreqs);
  }

  emxInit_boolean_T(&st, &x, 2, &y_emlrtRTEI, true);

  /*  all possible frequencies */
  /* Set max frequency to nyquist if only lower bound specified */
  /* Get just the frequencies for the given frequency range */
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, x, i, &t_emlrtRTEI);
  nz = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < nz; i++) {
    x->data[i] = (sfreqs->data[i] >= frequency_range[0]);
  }

  emxInit_boolean_T(sp, &r, 2, &w_emlrtRTEI, true);
  i = r->size[0] * r->size[1];
  r->size[0] = 1;
  r->size[1] = sfreqs->size[1];
  emxEnsureCapacity_boolean_T(sp, r, i, &u_emlrtRTEI);
  nz = sfreqs->size[0] * sfreqs->size[1];
  for (i = 0; i < nz; i++) {
    r->data[i] = (sfreqs->data[i] <= frequency_range[1]);
  }

  emxInit_boolean_T(sp, &freq_inds, 2, &v_emlrtRTEI, true);
  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])x->size, *(int32_T (*)[2])r->size,
    &emlrtECI, sp);
  i = freq_inds->size[0] * freq_inds->size[1];
  freq_inds->size[0] = 1;
  freq_inds->size[1] = x->size[1];
  emxEnsureCapacity_boolean_T(sp, freq_inds, i, &v_emlrtRTEI);
  nz = x->size[0] * x->size[1];
  for (i = 0; i < nz; i++) {
    freq_inds->data[i] = (x->data[i] && r->data[i]);
  }

  nz = x->size[1] - 1;
  trueCount = 0;
  for (b_i = 0; b_i <= nz; b_i++) {
    if (x->data[b_i] && r->data[b_i]) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (b_i = 0; b_i <= nz; b_i++) {
    if (x->data[b_i] && r->data[b_i]) {
      i = b_i + 1;
      if ((i < 1) || (i > sfreqs->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i, 1, sfreqs->size[1], &d_emlrtBCI, sp);
      }

      sfreqs->data[partialTrueCount] = sfreqs->data[i - 1];
      partialTrueCount++;
    }
  }

  i = sfreqs->size[0] * sfreqs->size[1];
  sfreqs->size[0] = 1;
  sfreqs->size[1] = trueCount;
  emxEnsureCapacity_real_T(sp, sfreqs, i, &w_emlrtRTEI);

  /* Compute the times of the middle of each spectrum */
  b = muDoubleScalarRound(winsize_samples / 2.0);
  i = stimes->size[0] * stimes->size[1];
  stimes->size[0] = 1;
  stimes->size[1] = window_start->size[1];
  emxEnsureCapacity_real_T(sp, stimes, i, &x_emlrtRTEI);
  nz = window_start->size[0] * window_start->size[1];
  for (i = 0; i < nz; i++) {
    stimes->data[i] = (window_start->data[i] + b) / Fs;
  }

  /* Preallocate spectrogram and slice data for efficient parallel computing */
  st.site = &d_emlrtRSI;
  i = x->size[0] * x->size[1];
  trueCount = x->size[0] * x->size[1];
  x->size[0] = 1;
  emxEnsureCapacity_boolean_T(&st, x, trueCount, &y_emlrtRTEI);
  nz = i - 1;
  for (i = 0; i <= nz; i++) {
    x->data[i] = (x->data[i] && r->data[i]);
  }

  emxFree_boolean_T(&r);
  b_st.site = &y_emlrtRSI;
  c_st.site = &ab_emlrtRSI;
  nz = combineVectorElements(&c_st, x);
  emxFree_boolean_T(&x);
  if (nz < 0) {
    emlrtNonNegativeCheckR2012b(nz, &b_emlrtDCI, sp);
  }

  i = mt_spectrogram->size[0] * mt_spectrogram->size[1];
  mt_spectrogram->size[0] = nz;
  mt_spectrogram->size[1] = window_start->size[1];
  emxEnsureCapacity_real32_T(sp, mt_spectrogram, i, &ab_emlrtRTEI);
  nz *= window_start->size[1];
  for (i = 0; i < nz; i++) {
    mt_spectrogram->data[i] = 0.0F;
  }

  nz = window_start->size[1] - 1;
  emlrtEnterParallelRegion(sp, omp_in_parallel());
  emlrtPushJmpBuf(sp, &emlrtJBStack);

#pragma omp parallel \
 num_threads(emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs())) \
 private(r1,fft_range,r2,data_segment,b_x,r3,b_data_segment,fft_data,b_fft_range,magnitude,mt_spectrum,b_emlrtJBEnviron,h_st,i1,b_window_start,b_data,loop_ub,d,end,y,exitg1,c_window_start,c_data_segment) \
 firstprivate(d_st,e_st,f_st,g_st,emlrtHadParallelError,d_window_start)

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
      g_st.prev = &f_st;
      g_st.tls = f_st.tls;
      h_st.prev = &g_st;
      h_st.tls = g_st.tls;
      emxInit_real32_T(&d_st, &r1, 2, &w_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &fft_range, 2, &kb_emlrtRTEI, true);
      emxInit_int32_T(&d_st, &r2, 2, &pb_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &data_segment, 2, &ib_emlrtRTEI, true);
      emxInit_boolean_T(&d_st, &b_x, 2, &fb_emlrtRTEI, true);
      emxInit_real_T(&d_st, &r3, 2, &w_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &b_data_segment, 2, &eb_emlrtRTEI, true);
      emxInit_creal32_T(&d_st, &fft_data, 2, &ob_emlrtRTEI, true);
      emxInit_creal32_T(&d_st, &b_fft_range, 2, &jb_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &magnitude, 2, &nb_emlrtRTEI, true);
      emxInit_real32_T(&d_st, &mt_spectrum, 1, &mb_emlrtRTEI, true);
    } else {
      emlrtHadParallelError = true;
    }

#pragma omp for nowait

    for (n = 0; n <= nz; n++) {
      if (emlrtHadParallelError)
        continue;
      if (setjmp(b_emlrtJBEnviron) == 0) {
        /*  window_idxs = window_start' + (0:winsize_samples-1); */
        /*  data_segments = data(window_idxs); */
        /*  COMPUTE THE MULTITAPER SPECTROGRAM */
        /*  */
        /*      STEP 1: Compute DPSS tapers based on desired spectral properties */
        /*      STEP 2: Multiply the data segment by the DPSS Tapers */
        /*      STEP 3: Compute the spectrum for each tapered segment */
        /*      STEP 4: Take the mean of the tapered spectra */
        /* Loop in parallel over all of the windows */
        /* Grab the data for the given window */
        if (muDoubleScalarIsNaN(winsize_samples - 1.0)) {
          i1 = r3->size[0] * r3->size[1];
          r3->size[0] = 1;
          r3->size[1] = 1;
          emxEnsureCapacity_real_T(&d_st, r3, i1, &db_emlrtRTEI);
          r3->data[0] = rtNaN;
        } else if (winsize_samples - 1.0 < 0.0) {
          r3->size[0] = 1;
          r3->size[1] = 0;
        } else if (muDoubleScalarIsInf(winsize_samples - 1.0) && (0.0 ==
                    winsize_samples - 1.0)) {
          i1 = r3->size[0] * r3->size[1];
          r3->size[0] = 1;
          r3->size[1] = 1;
          emxEnsureCapacity_real_T(&d_st, r3, i1, &db_emlrtRTEI);
          r3->data[0] = rtNaN;
        } else {
          i1 = r3->size[0] * r3->size[1];
          r3->size[0] = 1;
          loop_ub = (int32_T)muDoubleScalarFloor(winsize_samples - 1.0);
          r3->size[1] = loop_ub + 1;
          emxEnsureCapacity_real_T(&d_st, r3, i1, &db_emlrtRTEI);
          for (i1 = 0; i1 <= loop_ub; i1++) {
            r3->data[i1] = i1;
          }
        }

        e_st.site = &e_emlrtRSI;
        indexShapeCheck(&e_st, *(int32_T (*)[2])data->size, *(int32_T (*)[2])
                        r3->size);
        i1 = b_data_segment->size[0] * b_data_segment->size[1];
        b_data_segment->size[0] = 1;
        b_data_segment->size[1] = r3->size[1];
        emxEnsureCapacity_real32_T(&d_st, b_data_segment, i1, &eb_emlrtRTEI);
        i1 = n + 1;
        if ((i1 < 1) || (i1 > window_start->size[1])) {
          emlrtDynamicBoundsCheckR2012b(i1, 1, window_start->size[1],
            &c_emlrtBCI, &d_st);
        }

        b_window_start = window_start->data[i1 - 1];
        b_data = data->size[0] * data->size[1];
        loop_ub = r3->size[0] * r3->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          d = b_window_start + r3->data[i1];
          if (d != (int32_T)muDoubleScalarFloor(d)) {
            emlrtIntegerCheckR2012b(d, &emlrtDCI, &d_st);
          }

          end = (int32_T)d;
          if ((end < 1) || (end > b_data)) {
            emlrtDynamicBoundsCheckR2012b(end, 1, b_data, &b_emlrtBCI, &d_st);
          }

          b_data_segment->data[i1] = data->data[end - 1];
        }

        /* Skip empty segments */
        e_st.site = &f_emlrtRSI;
        i1 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = b_data_segment->size[1];
        emxEnsureCapacity_boolean_T(&e_st, b_x, i1, &fb_emlrtRTEI);
        loop_ub = b_data_segment->size[0] * b_data_segment->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          b_x->data[i1] = (b_data_segment->data[i1] == 0.0F);
        }

        f_st.site = &eb_emlrtRSI;
        y = true;
        g_st.site = &fb_emlrtRSI;
        if ((1 <= b_x->size[1]) && (b_x->size[1] > 2147483646)) {
          h_st.site = &q_emlrtRSI;
          check_forloop_overflow_error(&h_st);
        }

        b_data = 1;
        exitg1 = false;
        while ((!exitg1) && (b_data <= b_x->size[1])) {
          if (!b_x->data[b_data - 1]) {
            y = false;
            exitg1 = true;
          } else {
            b_data++;
          }
        }

        if (!y) {
          i1 = b_x->size[0] * b_x->size[1];
          b_x->size[0] = 1;
          b_x->size[1] = b_data_segment->size[1];
          emxEnsureCapacity_boolean_T(&d_st, b_x, i1, &gb_emlrtRTEI);
          loop_ub = b_data_segment->size[0] * b_data_segment->size[1];
          for (i1 = 0; i1 < loop_ub; i1++) {
            b_x->data[i1] = muSingleScalarIsNaN(b_data_segment->data[i1]);
          }

          e_st.site = &g_emlrtRSI;
          f_st.site = &gb_emlrtRSI;
          y = false;
          g_st.site = &fb_emlrtRSI;
          if ((1 <= b_x->size[1]) && (b_x->size[1] > 2147483646)) {
            h_st.site = &q_emlrtRSI;
            check_forloop_overflow_error(&h_st);
          }

          b_data = 1;
          exitg1 = false;
          while ((!exitg1) && (b_data <= b_x->size[1])) {
            if (!b_x->data[b_data - 1]) {
              b_data++;
            } else {
              y = true;
              exitg1 = true;
            }
          }

          if (y) {
            i1 = n + 1;
            if ((i1 < 1) || (i1 > mt_spectrogram->size[1])) {
              emlrtDynamicBoundsCheckR2012b(i1, 1, mt_spectrogram->size[1],
                &e_emlrtBCI, &d_st);
            }

            loop_ub = mt_spectrogram->size[0];
            for (i1 = 0; i1 < loop_ub; i1++) {
              mt_spectrogram->data[i1 + mt_spectrogram->size[0] * n] = rtNaNF;
            }
          } else {
            /* Option to detrend_opt data to remove low frequency DC component */
            if (detrend_opt == 1.0) {
              i1 = data_segment->size[0] * data_segment->size[1];
              data_segment->size[0] = 1;
              data_segment->size[1] = b_data_segment->size[1];
              emxEnsureCapacity_real32_T(&d_st, data_segment, i1, &hb_emlrtRTEI);
              loop_ub = b_data_segment->size[0] * b_data_segment->size[1] - 1;
              for (i1 = 0; i1 <= loop_ub; i1++) {
                data_segment->data[i1] = b_data_segment->data[i1];
              }

              e_st.site = &h_emlrtRSI;
              detrend(&e_st, data_segment, b_data_segment);
            } else {
              if (detrend_opt == 2.0) {
                i1 = data_segment->size[0] * data_segment->size[1];
                data_segment->size[0] = 1;
                data_segment->size[1] = b_data_segment->size[1];
                emxEnsureCapacity_real32_T(&d_st, data_segment, i1,
                  &ib_emlrtRTEI);
                loop_ub = b_data_segment->size[0] * b_data_segment->size[1] - 1;
                for (i1 = 0; i1 <= loop_ub; i1++) {
                  data_segment->data[i1] = b_data_segment->data[i1];
                }

                e_st.site = &i_emlrtRSI;
                b_detrend(&e_st, data_segment, b_data_segment);
              }
            }

            /* Multiply the data by the tapers (STEP 2) */
            c_window_start[0] = b_data_segment->size[1];
            c_data_segment = *b_data_segment;
            d_window_start[0] = c_window_start[0];
            c_data_segment.size = &d_window_start[0];
            c_data_segment.numDimensions = 1;
            e_st.site = &j_emlrtRSI;
            repmat(&e_st, &c_data_segment, num_tapers, magnitude);
            emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])magnitude->size, *(int32_T
              (*)[2])DPSS_tapers->size, &b_emlrtECI, &d_st);
            loop_ub = magnitude->size[0] * magnitude->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              magnitude->data[i1] *= (real32_T)DPSS_tapers->data[i1];
            }

            /* Compute the FFT (STEP 3) */
            e_st.site = &k_emlrtRSI;
            fft(&e_st, magnitude, nfft, fft_data);
            end = freq_inds->size[1] - 1;
            b_data = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (freq_inds->data[loop_ub]) {
                b_data++;
              }
            }

            i1 = r2->size[0] * r2->size[1];
            r2->size[0] = 1;
            r2->size[1] = b_data;
            emxEnsureCapacity_int32_T(&d_st, r2, i1, &w_emlrtRTEI);
            b_data = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (freq_inds->data[loop_ub]) {
                r2->data[b_data] = loop_ub + 1;
                b_data++;
              }
            }

            loop_ub = fft_data->size[1];
            i1 = b_fft_range->size[0] * b_fft_range->size[1];
            b_fft_range->size[0] = r2->size[1];
            b_fft_range->size[1] = fft_data->size[1];
            emxEnsureCapacity_creal32_T(&d_st, b_fft_range, i1, &jb_emlrtRTEI);
            for (i1 = 0; i1 < loop_ub; i1++) {
              b_data = r2->size[1];
              for (end = 0; end < b_data; end++) {
                if ((r2->data[end] < 1) || (r2->data[end] > fft_data->size[0]))
                {
                  emlrtDynamicBoundsCheckR2012b(r2->data[end], 1, fft_data->
                    size[0], &emlrtBCI, &d_st);
                }

                b_fft_range->data[end + b_fft_range->size[0] * i1].re =
                  fft_data->data[(r2->data[end] + fft_data->size[0] * i1) - 1].
                  re;
                if ((r2->data[end] < 1) || (r2->data[end] > fft_data->size[0]))
                {
                  emlrtDynamicBoundsCheckR2012b(r2->data[end], 1, fft_data->
                    size[0], &emlrtBCI, &d_st);
                }

                b_fft_range->data[end + b_fft_range->size[0] * i1].im =
                  fft_data->data[(r2->data[end] + fft_data->size[0] * i1) - 1].
                  im;
              }
            }

            /* Take the FFT magnitude (STEP 4) */
            i1 = fft_range->size[0] * fft_range->size[1];
            fft_range->size[0] = b_fft_range->size[0];
            fft_range->size[1] = b_fft_range->size[1];
            emxEnsureCapacity_real32_T(&d_st, fft_range, i1, &kb_emlrtRTEI);
            loop_ub = b_fft_range->size[0] * b_fft_range->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              fft_range->data[i1] = b_fft_range->data[i1].im;
            }

            e_st.site = &l_emlrtRSI;
            power(&e_st, fft_range, magnitude);
            i1 = fft_range->size[0] * fft_range->size[1];
            fft_range->size[0] = b_fft_range->size[0];
            fft_range->size[1] = b_fft_range->size[1];
            emxEnsureCapacity_real32_T(&d_st, fft_range, i1, &lb_emlrtRTEI);
            loop_ub = b_fft_range->size[0] * b_fft_range->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              fft_range->data[i1] = b_fft_range->data[i1].re;
            }

            e_st.site = &l_emlrtRSI;
            power(&e_st, fft_range, r1);
            emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])magnitude->size, *(int32_T
              (*)[2])r1->size, &c_emlrtECI, &d_st);
            loop_ub = magnitude->size[0] * magnitude->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              magnitude->data[i1] += r1->data[i1];
            }

            e_st.site = &m_emlrtRSI;
            f_st.site = &we_emlrtRSI;
            g_st.site = &ab_emlrtRSI;
            h_st.site = &xe_emlrtRSI;
            blockedSummation(&h_st, magnitude, magnitude->size[1], mt_spectrum);

            /* Add the spectrum to the spectrogram */
            i1 = n + 1;
            if ((i1 < 1) || (i1 > mt_spectrogram->size[1])) {
              emlrtDynamicBoundsCheckR2012b(i1, 1, mt_spectrogram->size[1],
                &f_emlrtBCI, &d_st);
            }

            c_window_start[0] = mt_spectrogram->size[0];
            emlrtSubAssignSizeCheckR2012b(&c_window_start[0], 1,
              &mt_spectrum->size[0], 1, &d_emlrtECI, &d_st);
            loop_ub = mt_spectrum->size[0];
            for (i1 = 0; i1 < loop_ub; i1++) {
              mt_spectrogram->data[i1 + mt_spectrogram->size[0] * n] =
                mt_spectrum->data[i1];
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
      emlrtHeapReferenceStackLeaveScope(&d_st, 11);
      emxFree_real32_T(&mt_spectrum);
      emxFree_real32_T(&magnitude);
      emxFree_creal32_T(&b_fft_range);
      emxFree_creal32_T(&fft_data);
      emxFree_real32_T(&b_data_segment);
      emxFree_real_T(&r3);
      emxFree_boolean_T(&b_x);
      emxFree_real32_T(&data_segment);
      emxFree_int32_T(&r2);
      emxFree_real32_T(&fft_range);
      emxFree_real32_T(&r1);
    }
  }

  emlrtPopJmpBuf(sp, &emlrtJBStack);
  emlrtExitParallelRegion(sp, omp_in_parallel());
  emxFree_boolean_T(&freq_inds);
  emxFree_real_T(&window_start);
  emxInit_real32_T(sp, &b_mt_spectrogram, 2, &bb_emlrtRTEI, true);

  /* Compute the mean FFT magnitude (STEP 4) */
  b = Fs * Fs;
  i = b_mt_spectrogram->size[0] * b_mt_spectrogram->size[1];
  b_mt_spectrogram->size[0] = mt_spectrogram->size[1];
  b_mt_spectrogram->size[1] = mt_spectrogram->size[0];
  emxEnsureCapacity_real32_T(sp, b_mt_spectrogram, i, &bb_emlrtRTEI);
  nz = mt_spectrogram->size[0];
  for (i = 0; i < nz; i++) {
    b_i = mt_spectrogram->size[1];
    for (trueCount = 0; trueCount < b_i; trueCount++) {
      b_mt_spectrogram->data[trueCount + b_mt_spectrogram->size[0] * i] =
        mt_spectrogram->data[i + mt_spectrogram->size[0] * trueCount] /
        (real32_T)b / (real32_T)num_tapers;
    }
  }

  i = mt_spectrogram->size[0] * mt_spectrogram->size[1];
  mt_spectrogram->size[0] = b_mt_spectrogram->size[0];
  mt_spectrogram->size[1] = b_mt_spectrogram->size[1];
  emxEnsureCapacity_real32_T(sp, mt_spectrogram, i, &cb_emlrtRTEI);
  nz = b_mt_spectrogram->size[0] * b_mt_spectrogram->size[1];
  for (i = 0; i < nz; i++) {
    mt_spectrogram->data[i] = b_mt_spectrogram->data[i];
  }

  emxFree_real32_T(&b_mt_spectrogram);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (multitaper_spectrogram_coder_mex.c) */
