// Rust port of the multitaper spectrogram core loop.
//
// Mirrors `calc_mts_batch` in
// /Users/Mike/code/labcode_main/multitaper/python/multitaper_spectrogram_python.py
// and the one-sided PSD scaling applied by `multitaper_spectrogram`.
//
// This library crate exposes a `compute_spectrogram` function usable from both
// the CLI binary (`src/main.rs`) and the optional PyO3 bindings (behind the
// `python` feature, enabled by default when building the cdylib).

// Allow stylistic lints across the crate — index-based loops are clearer
// than iterator chains for numerical code with multiple parallel arrays, and
// the manual abs/multiple-of patterns mirror the reference implementations.
#![allow(
    clippy::manual_abs_diff,
    clippy::manual_is_multiple_of,
    clippy::manual_memcpy,
    clippy::needless_range_loop,
)]

use ndarray::{s, Array1, Array2, ArrayView1, ArrayView2, Axis};
use rayon::prelude::*;
use realfft::RealFftPlanner;

pub mod dpss;
pub use dpss::{dpss, DpssError};

#[cfg(feature = "python")]
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
#[cfg(feature = "python")]
use pyo3::prelude::*;

// -------- Core types --------

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum DetrendMode {
    Linear,
    Constant,
    Off,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Weighting {
    Unity,
    Eigen,
}

impl DetrendMode {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "linear" => Ok(DetrendMode::Linear),
            "constant" | "const" => Ok(DetrendMode::Constant),
            "off" | "none" | "false" => Ok(DetrendMode::Off),
            other => Err(format!("Unknown detrend_opt: {}", other)),
        }
    }
}

impl Weighting {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "unity" => Ok(Weighting::Unity),
            "eigen" => Ok(Weighting::Eigen),
            other => Err(format!("Unsupported weighting: {}", other)),
        }
    }
}

// -------- Detrend helpers --------

/// In-place linear detrend of a 1-D slice (matches scipy.signal.detrend type='linear').
pub fn detrend_linear(row: &mut [f64]) {
    let n = row.len();
    if n < 2 {
        return;
    }
    let nf = n as f64;
    let x_mean = (nf - 1.0) * 0.5;
    let mut y_mean = 0.0;
    for &v in row.iter() {
        y_mean += v;
    }
    y_mean /= nf;

    let mut num = 0.0;
    let mut den = 0.0;
    for (i, &v) in row.iter().enumerate() {
        let dx = i as f64 - x_mean;
        num += dx * (v - y_mean);
        den += dx * dx;
    }
    let slope = if den != 0.0 { num / den } else { 0.0 };
    let intercept = y_mean - slope * x_mean;
    for (i, v) in row.iter_mut().enumerate() {
        *v -= intercept + slope * (i as f64);
    }
}

pub fn detrend_constant(row: &mut [f64]) {
    let n = row.len();
    if n == 0 {
        return;
    }
    let mean: f64 = row.iter().sum::<f64>() / (n as f64);
    for v in row.iter_mut() {
        *v -= mean;
    }
}

// -------- DPSS --------
//
// Native Rust DPSS solver in `dpss.rs` (pure-Rust port of
// `scipy.signal.windows.dpss(N, NW, K, return_ratios=True)`, validated
// against MATLAB R2025a `dpss(N, NW, K)` to <=1e-9). `compute_spectrogram`
// still accepts pre-computed `tapers` as input (shape `(K, winsize)`) for
// callers that already have them; new code can call `dpss::dpss(...)`.

/// Parameters controlling the compute.
#[derive(Clone, Debug)]
pub struct SpectrogramParams {
    pub fs: f64,
    /// Inclusive (min, max) frequency range.
    pub frequency_range: (f64, f64),
    /// Window (size_s, step_s) in seconds.
    pub window_params: (f64, f64),
    pub nfft: usize,
    pub detrend: DetrendMode,
    pub weighting: Weighting,
}

/// Outputs matching the Python function.
pub struct SpectrogramOutput {
    /// (nfreq_out, num_windows) PSD-scaled spectrogram.
    pub mt_spectrogram: Array2<f64>,
    /// (num_windows,) stimes in seconds.
    pub stimes: Array1<f64>,
    /// (nfreq_out,) sfreqs in Hz.
    pub sfreqs: Array1<f64>,
}

/// Compute the multitaper spectrogram given data, tapers, and parameters.
///
/// `data`    — 1-D float64 time series.
/// `tapers`  — shape (K, winsize), matching `scipy.signal.windows.dpss`.
/// `eigen`   — eigenvalues (length K), required when `weighting == Eigen`.
pub fn compute_spectrogram(
    data: ArrayView1<f64>,
    tapers: ArrayView2<f64>,
    eigen: Option<ArrayView1<f64>>,
    params: &SpectrogramParams,
) -> Result<SpectrogramOutput, String> {
    let fs = params.fs;
    let (fmin, fmax) = params.frequency_range;
    let winsize_samples = (params.window_params.0 * fs).round() as usize;
    let winstep_samples = (params.window_params.1 * fs).round() as usize;
    let nfft = params.nfft;

    if tapers.shape()[1] != winsize_samples {
        return Err(format!(
            "Taper length {} does not match window size {}",
            tapers.shape()[1],
            winsize_samples
        ));
    }
    let num_tapers = tapers.shape()[0];

    let wt: Array1<f64> = match params.weighting {
        Weighting::Unity => Array1::from_elem(num_tapers, 1.0 / num_tapers as f64),
        Weighting::Eigen => {
            let eigen =
                eigen.ok_or_else(|| "eigen is required when weighting='eigen'".to_string())?;
            if eigen.len() != num_tapers {
                return Err("eigen length mismatch with number of tapers".into());
            }
            eigen.mapv(|v| v / num_tapers as f64)
        }
    };

    let len_data = data.len();
    if len_data < winsize_samples {
        return Err("data shorter than window".into());
    }
    let mut window_starts: Vec<usize> = Vec::new();
    let mut start = 0usize;
    while start + winsize_samples <= len_data {
        window_starts.push(start);
        start += winstep_samples;
    }
    let num_windows = window_starts.len();

    // Frequency bins (one-sided, rfftfreq-style)
    let nfreq_full = nfft / 2 + 1;
    let df = fs / nfft as f64;
    let sfreqs_full: Vec<f64> = (0..nfreq_full).map(|k| k as f64 * df).collect();
    let freq_idx: Vec<usize> = sfreqs_full
        .iter()
        .enumerate()
        .filter_map(|(i, &f)| {
            if f >= fmin && f <= fmax {
                Some(i)
            } else {
                None
            }
        })
        .collect();
    let nfreq_out = freq_idx.len();
    let sfreqs_out: Array1<f64> = Array1::from_iter(freq_idx.iter().map(|&i| sfreqs_full[i]));

    // stimes — Python: (window_start + round(winsize/2)) / fs
    let half_win = (winsize_samples as f64 / 2.0).round();
    let stimes: Array1<f64> =
        Array1::from_iter(window_starts.iter().map(|&ws| (ws as f64 + half_win) / fs));

    // Transpose tapers to (winsize, K) to match Python's dpss_tapers_T layout.
    let tapers_t: Array2<f64> = tapers.t().to_owned();

    // Output (nfreq_out, num_windows) with PSD scaling applied directly.
    let mut output: Array2<f64> = Array2::<f64>::zeros((nfreq_out, num_windows));

    // PSD scale per output bin: /fs, x2 for non-DC/non-Nyquist.
    let nyquist_k = if nfft % 2 == 0 { Some(nfft / 2) } else { None };
    let bin_scale: Vec<f64> = freq_idx
        .iter()
        .map(|&k| {
            let is_dc = k == 0;
            let is_ny = Some(k) == nyquist_k;
            let m = if is_dc || is_ny { 1.0 } else { 2.0 };
            m / fs
        })
        .collect();

    let detrend_mode = params.detrend;
    let batch_size: usize = 1024;

    for chunk_start in (0..num_windows).step_by(batch_size) {
        let chunk_end = (chunk_start + batch_size).min(num_windows);
        let chunk_len = chunk_end - chunk_start;

        let tapers_t_view = tapers_t.view();
        let mut chunk_out: Array2<f64> = Array2::<f64>::zeros((chunk_len, nfreq_out));

        chunk_out
            .axis_iter_mut(Axis(0))
            .into_par_iter()
            .enumerate()
            .for_each_init(
                || {
                    let mut planner = RealFftPlanner::<f64>::new();
                    let f = planner.plan_fft_forward(nfft);
                    let in_buf = f.make_input_vec();
                    let out_buf = f.make_output_vec();
                    let scratch = f.make_scratch_vec();
                    (f, in_buf, out_buf, scratch)
                },
                |(f, in_buf, out_buf, scratch), (i_in_chunk, mut out_row)| {
                    let w = chunk_start + i_in_chunk;
                    let start_idx = window_starts[w];
                    let seg_slice = data.slice(s![start_idx..start_idx + winsize_samples]);
                    let mut seg: Vec<f64> = seg_slice.to_vec();

                    let all_zero = seg.iter().all(|&v| v == 0.0);
                    if all_zero {
                        return;
                    }

                    match detrend_mode {
                        DetrendMode::Linear => detrend_linear(&mut seg),
                        DetrendMode::Constant => detrend_constant(&mut seg),
                        DetrendMode::Off => {}
                    }

                    let mut power_accum: Vec<f64> = vec![0.0; nfreq_full];
                    for k in 0..num_tapers {
                        for t in 0..winsize_samples {
                            in_buf[t] = seg[t] * tapers_t_view[[t, k]];
                        }
                        for t in winsize_samples..nfft {
                            in_buf[t] = 0.0;
                        }
                        f.process_with_scratch(in_buf, out_buf, scratch)
                            .expect("rfft failed");

                        let wk = wt[k];
                        for j in 0..nfreq_full {
                            let c = out_buf[j];
                            let p = c.re * c.re + c.im * c.im;
                            power_accum[j] += wk * p;
                        }
                    }

                    for (oi, (&kk, &scale)) in
                        freq_idx.iter().zip(bin_scale.iter()).enumerate()
                    {
                        out_row[oi] = power_accum[kk] * scale;
                    }
                },
            );

        // Scatter chunk rows (chunk_len, nfreq_out) into output (nfreq_out, num_windows).
        for (i_in_chunk, row) in chunk_out.axis_iter(Axis(0)).enumerate() {
            let w = chunk_start + i_in_chunk;
            for (oi, v) in row.iter().enumerate() {
                output[[oi, w]] = *v;
            }
        }
    }

    Ok(SpectrogramOutput {
        mt_spectrogram: output,
        stimes,
        sfreqs: sfreqs_out,
    })
}

// -------- DPSS via linalg (Slepian sequences) --------
//
// We use the same approach scipy.signal.windows.dpss does: build a symmetric
// tridiagonal matrix and solve for the K largest eigenvectors. We call out to
// scipy's fortran-backed LAPACK via... no — we're in Rust. Reimplementing this
// here would mean pulling LAPACK bindings in, which is a heavy dependency we
// want to avoid. Instead the Python wrapper computes DPSS via scipy and passes
// tapers into `compute_spectrogram`. See `py_compute_spectrogram` below: it
// takes tapers as an input from the caller (numpy array), so Python retains
// the scipy DPSS call.

// -------- PyO3 bindings --------

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "compute_spectrogram", signature = (
    data,
    tapers,
    fs,
    frequency_range,
    window_params,
    nfft,
    detrend_opt = "linear".to_string(),
    weighting = "unity".to_string(),
    eigen = None,
))]
#[allow(clippy::too_many_arguments)]
fn py_compute_spectrogram<'py>(
    py: Python<'py>,
    data: PyReadonlyArray1<'py, f64>,
    tapers: numpy::PyReadonlyArray2<'py, f64>,
    fs: f64,
    frequency_range: (f64, f64),
    window_params: (f64, f64),
    nfft: usize,
    detrend_opt: String,
    weighting: String,
    eigen: Option<PyReadonlyArray1<'py, f64>>,
) -> PyResult<(
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
)> {
    let detrend = DetrendMode::parse(&detrend_opt)
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    let weighting_mode = Weighting::parse(&weighting)
        .map_err(pyo3::exceptions::PyValueError::new_err)?;

    let data_view = data.as_array();
    let tapers_view = tapers.as_array();
    let eigen_owned = eigen.as_ref().map(|e| e.as_array().to_owned());

    let params = SpectrogramParams {
        fs,
        frequency_range,
        window_params,
        nfft,
        detrend,
        weighting: weighting_mode,
    };

    let out = py.allow_threads(|| {
        compute_spectrogram(
            data_view,
            tapers_view,
            eigen_owned.as_ref().map(|e| e.view()),
            &params,
        )
    })
    .map_err(pyo3::exceptions::PyValueError::new_err)?;

    Ok((
        out.mt_spectrogram.into_pyarray_bound(py),
        out.stimes.into_pyarray_bound(py),
        out.sfreqs.into_pyarray_bound(py),
    ))
}

#[cfg(feature = "python")]
#[pymodule]
fn multitaper_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_compute_spectrogram, m)?)?;
    Ok(())
}
