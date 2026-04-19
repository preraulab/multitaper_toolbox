// CLI wrapper around `multitaper_rs::compute_spectrogram`.
//
// Reads data + tapers (+optional eigen) from .npy, computes the multitaper
// spectrogram, writes the result to .npy, and prints a TIMING line.

use std::time::Instant;

use clap::Parser;
use ndarray::{Array1, Array2};
use ndarray_npy::{ReadNpyExt, WriteNpyExt};

use multitaper_rs::{compute_spectrogram, DetrendMode, SpectrogramParams, Weighting};

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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let detrend = DetrendMode::parse(&cli.detrend_opt)?;
    let weighting = Weighting::parse(&cli.weighting)?;

    if cli.frequency_range.len() != 2 {
        return Err("--frequency-range expects two comma-separated values".into());
    }
    if cli.window_params.len() != 2 {
        return Err("--window-params expects two comma-separated values".into());
    }

    let params = SpectrogramParams {
        fs: cli.fs,
        frequency_range: (cli.frequency_range[0], cli.frequency_range[1]),
        window_params: (cli.window_params[0], cli.window_params[1]),
        nfft: cli.nfft,
        detrend,
        weighting,
    };

    let t0 = Instant::now();

    // Load data + tapers
    let data: Array1<f64> = Array1::<f64>::read_npy(std::fs::File::open(&cli.data)?)?;
    let tapers: Array2<f64> = Array2::<f64>::read_npy(std::fs::File::open(&cli.tapers)?)?;
    let eigen: Option<Array1<f64>> = match &cli.eigen {
        Some(path) => Some(Array1::<f64>::read_npy(std::fs::File::open(path)?)?),
        None => None,
    };

    let t_setup = t0.elapsed();

    let t_compute_start = Instant::now();
    let result = compute_spectrogram(
        data.view(),
        tapers.view(),
        eigen.as_ref().map(|e| e.view()),
        &params,
    )
    .map_err(|e| -> Box<dyn std::error::Error> { e.into() })?;
    let t_compute = t_compute_start.elapsed();

    let t_save_start = Instant::now();
    result
        .mt_spectrogram
        .write_npy(std::fs::File::create(&cli.output)?)?;
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
