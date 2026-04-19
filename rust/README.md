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
