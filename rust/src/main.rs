// Rust port of the multitaper spectrogram core loop.
//
// Mirrors `calc_mts_batch` in
// /Users/Mike/code/labcode_main/multitaper/python/multitaper_spectrogram_python.py
// and the one-sided PSD scaling applied by `multitaper_spectrogram`.
//
// DPSS tapers are NOT computed here -- they must be supplied via an .npy
// file produced by `scipy.signal.windows.dpss`. 'adapt' weighting is not
// implemented (follow-up).

use std::time::Instant;

use clap::Parser;
use ndarray::{s, Array1, Array2, Axis};
use ndarray_npy::{ReadNpyExt, WriteNpyExt};
use rayon::prelude::*;
use realfft::RealFftPlanner;

#[derive(Parser, Debug)]
#[command(author, version, about = "Multitaper spectrogram (Rust)")]
struct Cli {
    /// Path to .npy file containing 1-D float64 data.
    #[arg(long)]
    data: String,

    /// Path to .npy file containing (K, winsize) float64 DPSS tapers.
    #[arg(long)]
    tapers: String,

    /// Optional path to .npy file with float64 DPSS eigenvalues (length K).
    /// Required only when --weighting eigen.
    #[arg(long)]
    eigen: Option<String>,

    /// Sampling frequency (Hz).
    #[arg(long)]
    fs: f64,

    /// Frequency range as "min,max" (Hz).
    #[arg(long, value_delimiter = ',')]
    frequency_range: Vec<f64>,

    /// Window params as "size_s,step_s" (seconds).
    #[arg(long, value_delimiter = ',')]
    window_params: Vec<f64>,

    /// FFT length (nfft).
    #[arg(long)]
    nfft: usize,

    /// Detrend option: linear | constant | off.
    #[arg(long, default_value = "linear")]
    detrend_opt: String,

    /// Weighting: unity | eigen.
    #[arg(long, default_value = "unity")]
    weighting: String,

    /// Output .npy path (float64 mt_spectrogram, shape (nfreq, num_windows)).
    #[arg(long)]
    output: String,

    /// Print timing info to stderr.
    #[arg(long, default_value_t = false)]
    verbose: bool,
}

/// In-place linear detrend of a 1-D slice (matches scipy.signal.detrend type='linear').
fn detrend_linear(row: &mut [f64]) {
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

fn detrend_constant(row: &mut [f64]) {
    let n = row.len();
    if n == 0 {
        return;
    }
    let mean: f64 = row.iter().sum::<f64>() / (n as f64);
    for v in row.iter_mut() {
        *v -= mean;
    }
}

#[derive(Copy, Clone, PartialEq, Eq)]
enum DetrendMode {
    Linear,
    Constant,
    Off,
}

#[derive(Copy, Clone, PartialEq, Eq)]
enum Weighting {
    Unity,
    Eigen,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let detrend_mode = match cli.detrend_opt.as_str() {
        "linear" => DetrendMode::Linear,
        "constant" | "const" => DetrendMode::Constant,
        "off" | "none" | "false" => DetrendMode::Off,
        other => return Err(format!("Unknown detrend_opt: {}", other).into()),
    };
    let weighting = match cli.weighting.as_str() {
        "unity" => Weighting::Unity,
        "eigen" => Weighting::Eigen,
        other => return Err(format!("Unsupported weighting: {}", other).into()),
    };

    if cli.frequency_range.len() != 2 {
        return Err("--frequency-range expects two comma-separated values".into());
    }
    if cli.window_params.len() != 2 {
        return Err("--window-params expects two comma-separated values".into());
    }
    let frange = (cli.frequency_range[0], cli.frequency_range[1]);
    let fs = cli.fs;
    let winsize_samples = (cli.window_params[0] * fs).round() as usize;
    let winstep_samples = (cli.window_params[1] * fs).round() as usize;
    let nfft = cli.nfft;

    let t0 = Instant::now();

    // Load data + tapers
    let data: Array1<f64> = Array1::<f64>::read_npy(std::fs::File::open(&cli.data)?)?;
    let tapers: Array2<f64> = Array2::<f64>::read_npy(std::fs::File::open(&cli.tapers)?)?;
    if tapers.shape()[1] != winsize_samples {
        return Err(format!(
            "Taper length {} does not match window size {}",
            tapers.shape()[1],
            winsize_samples
        )
        .into());
    }
    let num_tapers = tapers.shape()[0];

    // Weights vector (length K); 'unity' -> 1/K, 'eigen' -> eigen / K.
    let wt: Array1<f64> = match weighting {
        Weighting::Unity => Array1::from_elem(num_tapers, 1.0 / num_tapers as f64),
        Weighting::Eigen => {
            let eigen_path = cli
                .eigen
                .as_ref()
                .ok_or("--eigen is required when --weighting eigen")?;
            let eigen: Array1<f64> =
                Array1::<f64>::read_npy(std::fs::File::open(eigen_path)?)?;
            if eigen.len() != num_tapers {
                return Err("eigen length mismatch with number of tapers".into());
            }
            &eigen / (num_tapers as f64)
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
    let freq_mask: Vec<bool> = sfreqs_full
        .iter()
        .map(|&f| f >= frange.0 && f <= frange.1)
        .collect();
    let freq_idx: Vec<usize> = freq_mask
        .iter()
        .enumerate()
        .filter_map(|(i, &m)| if m { Some(i) } else { None })
        .collect();
    let nfreq_out = freq_idx.len();

    // Transpose tapers to (winsize, K) to match Python's dpss_tapers_T layout.
    let tapers_t: Array2<f64> = tapers.t().to_owned();

    if cli.verbose {
        eprintln!(
            "[rust] num_windows={} winsize={} K={} nfft={} nfreq_out={}",
            num_windows, winsize_samples, num_tapers, nfft, nfreq_out
        );
    }

    let t_setup = t0.elapsed();
    let t_compute_start = Instant::now();

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

    let t_compute = t_compute_start.elapsed();

    let t_save_start = Instant::now();
    output.write_npy(std::fs::File::create(&cli.output)?)?;
    let t_save = t_save_start.elapsed();

    if cli.verbose {
        eprintln!(
            "[rust] setup={:?} compute={:?} save={:?} total={:?}",
            t_setup,
            t_compute,
            t_save,
            t0.elapsed()
        );
    }
    println!(
        "TIMING setup_s={:.6} compute_s={:.6} save_s={:.6} total_s={:.6}",
        t_setup.as_secs_f64(),
        t_compute.as_secs_f64(),
        t_save.as_secs_f64(),
        t0.elapsed().as_secs_f64()
    );

    Ok(())
}
