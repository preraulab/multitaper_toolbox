"""Long-signal speedup benchmark: Python vs Rust multitaper.

Runs `multitaper_spectrogram` on increasingly long signals at Fs=200 Hz and
reports median wall-clock time over 3 trials (one warm-up run is discarded)
for both the pure-Python and Rust backends.

Usage:
    python test_rust_speedup.py

Prints a compact markdown table at the end:

    | Duration | Python (s) | Rust (s) | Ratio |
"""
from __future__ import annotations

import gc
import os
import platform
import sys
import time
from statistics import median

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from multitaper_spectrogram_python import (  # noqa: E402
    multitaper_spectrogram,
    _HAS_RUST,
)


def make_signal(duration_s: float, fs: float, seed: int = 1) -> np.ndarray:
    """Same generator as test_rust_equivalence.py (fixed seed)."""
    rng = np.random.default_rng(seed)
    n = int(round(duration_s * fs))
    t = np.arange(n) / fs
    sig = (
        1.0 * np.sin(2 * np.pi * 2.0 * t)
        + 0.7 * np.sin(2 * np.pi * 7.5 * t)
        + 0.4 * np.sin(2 * np.pi * 18.0 * t)
        + 0.5 * rng.standard_normal(n)
    )
    return sig.astype(np.float64)


def run_one(data, fs, params, use_rust):
    return multitaper_spectrogram(
        data,
        fs,
        frequency_range=params["frequency_range"],
        time_bandwidth=params["time_bandwidth"],
        num_tapers=params["num_tapers"],
        window_params=params["window_params"],
        min_nfft=params.get("min_nfft", 0),
        detrend_opt=params.get("detrend_opt", "linear"),
        weighting=params.get("weighting", "unity"),
        plot_on=False,
        verbose=False,
        use_rust=use_rust,
    )


def time_trials(data, fs, params, use_rust, n_trials=3):
    # Warm-up (discarded).
    _ = run_one(data, fs, params, use_rust=use_rust)
    gc.collect()
    times = []
    for _ in range(n_trials):
        t0 = time.perf_counter()
        _ = run_one(data, fs, params, use_rust=use_rust)
        times.append(time.perf_counter() - t0)
        gc.collect()
    return median(times), times


def main():
    if not _HAS_RUST:
        print("ERROR: multitaper_rs not installed. Aborting.")
        sys.exit(1)

    print(f"Platform: {platform.machine()} / {platform.processor()}")
    print(f"Python:   {platform.python_version()} @ {sys.executable}")
    print(f"numpy:    {np.__version__}")

    fs = 200.0
    params = dict(
        frequency_range=[0, 25],
        time_bandwidth=5,
        num_tapers=9,
        window_params=[5, 1],
        detrend_opt="linear",
        weighting="unity",
    )

    # Durations in hours. 24h will be skipped if signal doesn't fit.
    duration_hours = [1, 4, 10, 24]

    # Very rough memory check: float64 samples plus scratch ~ 3x signal size.
    # 24h @ 200 Hz = 17.28M samples * 8 B = ~138 MB; trivially fits.
    # Bigger concern is the spectrogram output + internal buffers.

    rows = []
    for hrs in duration_hours:
        dur_s = hrs * 3600.0
        label = f"{hrs} h"
        print(f"\n=== {label} @ {fs} Hz ===")
        try:
            data = make_signal(dur_s, fs, seed=1)
        except MemoryError:
            print(f"    skipped — MemoryError allocating signal")
            rows.append((label, None, None))
            continue
        print(f"    samples={data.size:,}  ({data.nbytes/1e6:.1f} MB)")

        try:
            t_py, trials_py = time_trials(data, fs, params, use_rust=False)
            print(f"    python  median={t_py:.3f}s   trials={['%.3f'%x for x in trials_py]}")
        except MemoryError:
            print(f"    python  skipped — MemoryError")
            t_py = None

        try:
            t_rs, trials_rs = time_trials(data, fs, params, use_rust=True)
            print(f"    rust    median={t_rs:.3f}s   trials={['%.3f'%x for x in trials_rs]}")
        except MemoryError:
            print(f"    rust    skipped — MemoryError")
            t_rs = None

        rows.append((label, t_py, t_rs))

        del data
        gc.collect()

    # Compact markdown table.
    print("\n\n| Duration | Python (s) | Rust (s) | Ratio |")
    print("|----------|------------|----------|-------|")
    for label, t_py, t_rs in rows:
        if t_py is None or t_rs is None:
            py_s = "—" if t_py is None else f"{t_py:.2f}"
            rs_s = "—" if t_rs is None else f"{t_rs:.2f}"
            ratio = "—"
        else:
            py_s = f"{t_py:.2f}"
            rs_s = f"{t_rs:.2f}"
            ratio = f"{t_py / t_rs:.2f}x"
        print(f"| {label:>8} | {py_s:>10} | {rs_s:>8} | {ratio:>5} |")


if __name__ == "__main__":
    main()
