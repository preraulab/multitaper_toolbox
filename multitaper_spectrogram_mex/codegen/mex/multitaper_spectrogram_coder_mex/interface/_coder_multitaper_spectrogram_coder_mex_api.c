/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_multitaper_spectrogram_coder_mex_api.c
 *
 * Code generation for function '_coder_multitaper_spectrogram_coder_mex_api'
 *
 */

/* Include files */
#include "_coder_multitaper_spectrogram_coder_mex_api.h"
#include "multitaper_spectrogram_coder_mex.h"
#include "multitaper_spectrogram_coder_mex_data.h"
#include "multitaper_spectrogram_coder_mex_emxutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo wc_emlrtRTEI = { 1,/* lineNo */
  1,                                   /* colNo */
  "_coder_multitaper_spectrogram_coder_mex_api",/* fName */
  ""                                   /* pName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real32_T *y);
static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Fs, const
  char_T *identifier);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *frequency_range, const char_T *identifier))[2];
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *data, const
  char_T *identifier, emxArray_real32_T *y);
static const mxArray *emlrt_marshallOut(const emxArray_real32_T *u);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2];
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *DPSS_tapers,
  const char_T *identifier, emxArray_real_T *y);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real32_T *ret);
static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2];
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real32_T *y)
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u->data[0]);
  emlrtSetDimensions((mxArray *)m, u->size, 2);
  emlrtAssign(&y, m);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Fs, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(Fs), &thisId);
  emlrtDestroyArray(&Fs);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *frequency_range, const char_T *identifier))[2]
{
  real_T (*y)[2];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(frequency_range), &thisId);
  emlrtDestroyArray(&frequency_range);
  return y;
}
  static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *data, const
  char_T *identifier, emxArray_real32_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(data), &thisId, y);
  emlrtDestroyArray(&data);
}

static const mxArray *emlrt_marshallOut(const emxArray_real32_T *u)
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxSINGLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u->data[0]);
  emlrtSetDimensions((mxArray *)m, u->size, 2);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2]
{
  real_T (*y)[2];
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *DPSS_tapers, const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  h_emlrt_marshallIn(sp, emlrtAlias(DPSS_tapers), &thisId, y);
  emlrtDestroyArray(&DPSS_tapers);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real32_T *ret)
{
  static const int32_T dims[2] = { -1, -1 };

  const boolean_T bv[2] = { true, true };

  int32_T iv[2];
  int32_T i;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "single", false, 2U, dims, &bv[0],
    iv);
  ret->allocatedSize = iv[0] * iv[1];
  i = ret->size[0] * ret->size[1];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  emxEnsureCapacity_real32_T(sp, ret, i, (emlrtRTEInfo *)NULL);
  ret->data = (real32_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2]
{
  real_T (*ret)[2];
  static const int32_T dims[2] = { 1, 2 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[2])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { -1, -1 };

  const boolean_T bv[2] = { true, true };

  int32_T iv[2];
  int32_T i;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv[0],
    iv);
  ret->allocatedSize = iv[0] * iv[1];
  i = ret->size[0] * ret->size[1];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  emxEnsureCapacity_real_T(sp, ret, i, (emlrtRTEInfo *)NULL);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

void multitaper_spectrogram_coder_mex_api(const mxArray * const prhs[10],
  int32_T nlhs, const mxArray *plhs[3])
{
  emxArray_real32_T *data;
  emxArray_real_T *DPSS_tapers;
  emxArray_real32_T *mt_spectrogram;
  emxArray_real_T *stimes;
  emxArray_real_T *sfreqs;
  const mxArray *prhs_copy_idx_2;
  real_T Fs;
  real_T (*frequency_range)[2];
  real_T time_bandwidth;
  real_T num_tapers;
  real_T winsize_samples;
  real_T winstep_samples;
  real_T min_NFFT;
  real_T detrend_opt;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real32_T(&st, &data, 2, &wc_emlrtRTEI, true);
  emxInit_real_T(&st, &DPSS_tapers, 2, &wc_emlrtRTEI, true);
  emxInit_real32_T(&st, &mt_spectrogram, 2, &wc_emlrtRTEI, true);
  emxInit_real_T(&st, &stimes, 2, &wc_emlrtRTEI, true);
  emxInit_real_T(&st, &sfreqs, 2, &wc_emlrtRTEI, true);
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  data->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "data", data);
  Fs = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "Fs");
  frequency_range = e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2),
    "frequency_range");
  DPSS_tapers->canFreeData = false;
  g_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "DPSS_tapers", DPSS_tapers);
  time_bandwidth = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]),
    "time_bandwidth");
  num_tapers = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "num_tapers");
  winsize_samples = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]),
    "winsize_samples");
  winstep_samples = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]),
    "winstep_samples");
  min_NFFT = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "min_NFFT");
  detrend_opt = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "detrend_opt");

  /* Invoke the target function */
  multitaper_spectrogram_coder_mex(&st, data, Fs, *frequency_range, DPSS_tapers,
    time_bandwidth, num_tapers, winsize_samples, winstep_samples, min_NFFT,
    detrend_opt, mt_spectrogram, stimes, sfreqs);

  /* Marshall function outputs */
  mt_spectrogram->canFreeData = false;
  plhs[0] = emlrt_marshallOut(mt_spectrogram);
  emxFree_real32_T(&mt_spectrogram);
  emxFree_real_T(&DPSS_tapers);
  emxFree_real32_T(&data);
  if (nlhs > 1) {
    stimes->canFreeData = false;
    plhs[1] = b_emlrt_marshallOut(stimes);
  }

  emxFree_real_T(&stimes);
  if (nlhs > 2) {
    sfreqs->canFreeData = false;
    plhs[2] = b_emlrt_marshallOut(sfreqs);
  }

  emxFree_real_T(&sfreqs);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_multitaper_spectrogram_coder_mex_api.c) */
