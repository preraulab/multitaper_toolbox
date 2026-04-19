"""Equivalence + timing tests for the Rust backend.

Usage:
    python test_rust_equivalence.py

Runs multitaper_spectrogram twice for each of several parameter sets — once
with use_rust=False (pure Python) and once with use_rust=True (Rust backend) —
and asserts that outputs agree to a small tolerance. Also runs a one-hour
timing comparison at the end.
"""
from __future__ import annotations

import sys
import time
import numpy as np

# Local import so the script runs from its own directory without installing.
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from multitaper_spectrogram_python import (  # noqa: E402
    multitaper_spectrogram,
    _HAS_RUST,
)


def make_signal(duration_s: float, fs: float, seed: int = 0) -> np.ndarray:
    """Reproducible test signal: coloured noise + a few sinusoids."""
    rng = np.random.default_rng(seed)
    n = int(round(duration_s * fs))
    t = np.arange(n) / fs
    sig = (
        1.0 * np.sin(2 * np.pi * 2.0 * t)      # 2 Hz
        + 0.7 * np.sin(2 * np.pi * 7.5 * t)    # 7.5 Hz
        + 0.4 * np.sin(2 * np.pi * 18.0 * t)   # 18 Hz
        + 0.5 * rng.standard_normal(n)         # white noise
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


def check_equivalence(data, fs, params, label, atol=1e-8, rtol=1e-6):
    print(f"\n=== Equivalence test: {label} ===")
    print(f"    params={params}")

    t0 = time.perf_counter()
    py_out = run_one(data, fs, params, use_rust=False)
    t_py = time.perf_counter() - t0

    t0 = time.perf_counter()
    rs_out = run_one(data, fs, params, use_rust=True)
    t_rs = time.perf_counter() - t0

    py_spec, py_stimes, py_sfreqs = py_out
    rs_spec, rs_stimes, rs_sfreqs = rs_out

    print(f"    shape py={py_spec.shape}  rs={rs_spec.shape}")
    print(f"    python time={t_py:.3f}s   rust time={t_rs:.3f}s")

    assert py_spec.shape == rs_spec.shape, (
        f"shape mismatch py={py_spec.shape} rs={rs_spec.shape}"
    )
    assert np.allclose(py_stimes, rs_stimes, atol=1e-10), "stimes differ"
    assert np.allclose(py_sfreqs, rs_sfreqs, atol=1e-10), "sfreqs differ"

    abs_diff = np.abs(py_spec - rs_spec)
    max_abs = float(abs_diff.max())
    ref = np.maximum(np.abs(py_spec), np.abs(rs_spec))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_diff = np.where(ref > 0, abs_diff / ref, 0.0)
    max_rel = float(rel_diff.max())

    print(f"    max abs diff = {max_abs:.3e}")
    print(f"    max rel diff = {max_rel:.3e}")

    ok = np.allclose(py_spec, rs_spec, atol=atol, rtol=rtol)
    print(f"    allclose(atol={atol}, rtol={rtol}) = {ok}")
    if not ok:
        print(f"    FAILED tolerance — reporting numbers anyway")
    return ok, max_abs, max_rel


def timing_bench(duration_s: float, fs: float, params: dict):
    print(f"\n=== Timing benchmark: {duration_s/60:.1f} minutes @ {fs} Hz ===")
    data = make_signal(duration_s, fs, seed=1)

    # Warm-up run (discard).
    _ = run_one(data, fs, params, use_rust=False)
    _ = run_one(data, fs, params, use_rust=True)

    t0 = time.perf_counter()
    _ = run_one(data, fs, params, use_rust=False)
    t_py = time.perf_counter() - t0

    t0 = time.perf_counter()
    _ = run_one(data, fs, params, use_rust=True)
    t_rs = time.perf_counter() - t0

    print(f"    python: {t_py*1000:.1f} ms")
    print(f"    rust:   {t_rs*1000:.1f} ms")
    print(f"    speedup: {t_py / t_rs:.2f}x")
    return t_py, t_rs


def main():
    if not _HAS_RUST:
        print("ERROR: multitaper_rs extension is not installed.")
        print("Build with: cd multitaper/rust && maturin develop --release")
        sys.exit(1)

    fs = 200.0
    # 5-minute signal for equivalence checks
    data = make_signal(5 * 60, fs, seed=0)

    param_sets = [
        (
            "default 5s/1s window, unity weighting",
            dict(
                frequency_range=[0, 25],
                time_bandwidth=5,
                num_tapers=9,
                window_params=[5, 1],
                detrend_opt="linear",
                weighting="unity",
            ),
        ),
        (
            "wider 10s/2s window, eigen weighting",
            dict(
                frequency_range=[0, 50],
                time_bandwidth=4,
                num_tapers=7,
                window_params=[10, 2],
                detrend_opt="constant",
                weighting="eigen",
            ),
        ),
        (
            "narrow 2s/0.5s window, detrend off",
            dict(
                frequency_range=[1, 40],
                time_bandwidth=3,
                num_tapers=5,
                window_params=[2, 0.5],
                detrend_opt="off",
                weighting="unity",
            ),
        ),
    ]

    results = []
    for label, params in param_sets:
        ok, ma, mr = check_equivalence(data, fs, params, label)
        results.append((label, ok, ma, mr))

    print("\n=== Equivalence summary ===")
    all_ok = True
    for label, ok, ma, mr in results:
        status = "PASS" if ok else "FAIL"
        print(f"    [{status}] {label}: max_abs={ma:.3e}  max_rel={mr:.3e}")
        all_ok = all_ok and ok

    # Timing benchmark: one hour @ 200 Hz.
    t_py, t_rs = timing_bench(
        duration_s=3600.0,
        fs=200.0,
        params=dict(
            frequency_range=[0, 25],
            time_bandwidth=5,
            num_tapers=9,
            window_params=[5, 1],
            detrend_opt="linear",
            weighting="unity",
        ),
    )

    print("\n=== Final ===")
    print(f"Equivalence: {'PASS' if all_ok else 'FAIL'}")
    print(f"Speedup (1h @ 200Hz): {t_py/t_rs:.2f}x")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
