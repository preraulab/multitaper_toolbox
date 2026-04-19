# multitaper_rs — Rust port of the multitaper spectrogram

This is a Rust implementation of the core multitaper spectrogram compute loop
from `../python/multitaper_spectrogram_python.py` (function `calc_mts_batch`,
plus the one-sided PSD scaling applied by `multitaper_spectrogram`).

## Scope

Implemented:
- Linear / constant / off detrending (per window).
- Taper multiply, zero-padded real FFT (pocketfft via `realfft`), weighted
  power average across tapers.
- One-sided PSD scaling (2x on non-DC/non-Nyquist bins, division by `fs`).
- Rayon parallelism over windows within chunks.
- `'unity'` and `'eigen'` taper weighting.

Not implemented:
- DPSS tapers. **Tapers must be pre-computed** (e.g. with
  `scipy.signal.windows.dpss`) and supplied via `.npy`. Reimplementing the
  Slepian eigenvalue solve in Rust is out of scope for this port; the FFT +
  power + weighting pipeline is what dominates runtime.
- `'adapt'` taper weighting (iterative, per-window) — follow-up.

## Build

```
cargo build --release
```

Binary is produced at `target/release/multitaper_rs`.

## Usage

```
multitaper_rs \
  --data data.npy --tapers tapers.npy \
  [--eigen eigen.npy] \
  --fs 200 --frequency-range 0,25 \
  --window-params 5,1 --nfft 1024 \
  --detrend-opt linear --weighting unity \
  --output output.npy
```

Output `.npy` has shape `(nfreq, num_windows)`, matching the Python
`mt_spectrogram` array that is returned after the one-sided PSD scaling block.

Timing info is printed to stdout as `TIMING setup_s=... compute_s=... save_s=... total_s=...`.

## Generating tapers

```python
from scipy.signal.windows import dpss
import numpy as np
tapers, eigen = dpss(winsize_samples, time_bandwidth, num_tapers, return_ratios=True)
np.save("tapers.npy", tapers)          # shape (K, winsize)
np.save("eigen.npy", eigen)            # only needed for --weighting eigen
```

## Python bindings

The same crate also builds as a Python extension (`multitaper_rs`) that
`../python/multitaper_spectrogram_python.py` will transparently use when
available. On machines without the extension installed the Python module
falls back to its pure-Python path automatically.

Build + install into the active Python environment:

```
cd multitaper/rust
pip install maturin   # once per env, if not already on PATH
maturin develop --release
```

`maturin develop --release` builds the cdylib with the `python` feature
enabled (see `pyproject.toml`) and installs it as an editable package into
the current interpreter. The crate is configured with `crate-type =
["cdylib", "rlib"]` in `Cargo.toml` so the same source powers both the
Python extension and the standalone CLI binary (`src/main.rs` links
against it as an `rlib`). The CLI binary is unaffected by the Python
build — `cargo build --release --bin multitaper_rs` still produces the
standalone tool with no Python dependencies.

### Using the backend from Python

```python
from multitaper_spectrogram_python import multitaper_spectrogram
# use_rust=None (default) -> auto-select Rust when installed
spec, stimes, sfreqs = multitaper_spectrogram(data, fs, use_rust=None)

# Force pure Python (A/B comparisons):
spec, stimes, sfreqs = multitaper_spectrogram(data, fs, use_rust=False)
```

Once per process a banner is printed to stderr:

```
multitaper: using Rust backend
multitaper: Rust backend unavailable, using pure Python
```

### Equivalence + timing test

```
python multitaper/python/test_rust_equivalence.py
```

Checks that Rust and pure-Python paths agree to `atol=1e-8, rtol=1e-6`
across three parameter sets, then benchmarks a 1-hour signal at 200 Hz.

### Limitations

- `weighting='adapt'` is not implemented in Rust; the Python wrapper
  falls back automatically for that case.
- DPSS tapers themselves are still computed by `scipy.signal.windows.dpss`
  in Python and passed into Rust. (Porting the Slepian eigensolver would
  require pulling in LAPACK.)
