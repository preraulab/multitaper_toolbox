//! Discrete Prolate Spheroidal Sequences (DPSS / Slepian) taper generation.
//!
//! Pure-Rust port of `scipy.signal.windows.dpss(N, NW, K, return_ratios=True)`.
//! Removes the MATLAB / scipy dependency for the standalone CLI.
//!
//! # Algorithm
//!
//! Slepian (1978) showed that the DPSS tapers are the eigenvectors of the
//! symmetric tridiagonal matrix `T` with
//!
//! ```text
//! diag(i)     = ((N - 1 - 2i) / 2)^2 * cos(2π · NW / N)         i = 0..N-1
//! off-diag(i) = (i + 1) · (N - 1 - i) / 2                       i = 0..N-2
//! ```
//!
//! corresponding to the K largest eigenvalues (Percival & Walden 1993, §8.3).
//! Crucially, the eigenvalues of `T` are *not* the concentration ratios — `T`
//! has the same eigenvectors as the prolate concentration matrix but
//! different eigenvalues. We compute the actual concentration ratios
//! separately via direct quadrature against the sinc kernel — matches
//! scipy's `return_ratios=True` branch.
//!
//! # Eigensolver
//!
//! Classical implicit-shift symmetric-tridiagonal QL (LAPACK DSTEQR /
//! Numerical Recipes `tqli`) with Wilkinson shifts. For DYNAM-O window
//! sizes (N ≤ 500) the O(N^2) work is negligible — DPSS is computed once
//! per (N, NW, K) tuple, dwarfed by the watershed cost.
//!
//! # Sign convention (scipy-compatible)
//!
//! - Even-index tapers (k = 0, 2, 4, …) are negated if their sum is negative,
//!   so the sum is non-negative.
//! - Odd-index tapers (k = 1, 3, 5, …) are negated if the first sample whose
//!   magnitude exceeds `max(1e-7, 1/N)` is negative.

use ndarray::{Array1, Array2};

#[derive(Debug, thiserror::Error)]
pub enum DpssError {
    #[error("invalid parameters: {0}")]
    InvalidParams(String),
    #[error("symmetric-tridiagonal eigensolver failed to converge after {0} iterations")]
    NoConvergence(usize),
}

/// Compute K DPSS tapers of length N with time-half-bandwidth NW.
///
/// Returns `(tapers, ratios)` where:
/// - `tapers` has shape `(K, N)`; rows are unit-L²-norm Slepian sequences.
/// - `ratios` has shape `(K,)`; entries are concentration ratios in `(0, 1]`,
///   monotonically decreasing.
pub fn dpss(n: usize, nw: f64, k: usize) -> Result<(Array2<f64>, Array1<f64>), DpssError> {
    if n < 2 {
        return Err(DpssError::InvalidParams(format!("N must be >= 2, got {}", n)));
    }
    if k == 0 || k > n {
        return Err(DpssError::InvalidParams(format!(
            "K must be in 1..=N (={}), got {}", n, k
        )));
    }
    if !(nw > 0.0 && nw < n as f64 / 2.0) {
        return Err(DpssError::InvalidParams(format!(
            "NW must be in (0, N/2), got NW={} for N={}", nw, n
        )));
    }

    // 1. Build the tridiagonal prolate matrix T (Slepian 1978).
    let nf = n as f64;
    let cw = (2.0 * std::f64::consts::PI * nw / nf).cos();
    let mut diag = vec![0.0_f64; n];
    let mut off = vec![0.0_f64; n - 1];
    for i in 0..n {
        let m = (nf - 1.0 - 2.0 * i as f64) / 2.0;
        diag[i] = m * m * cw;
    }
    for i in 0..n - 1 {
        let ip = i as f64 + 1.0;
        off[i] = ip * (nf - ip) / 2.0;
    }

    // 2. Symmetric-tridiagonal eigendecomposition. Returns eigenvalues sorted
    //    ascending and eigenvectors as columns of `vecs`.
    let (vals, vecs) = symmetric_tridiagonal_eigh(&diag, &off)?;

    // 3. Pick the K eigenpairs with the LARGEST eigenvalues — these correspond
    //    to the most-concentrated Slepian sequences.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| vals[b].partial_cmp(&vals[a]).unwrap());

    let mut tapers = Array2::<f64>::zeros((k, n));
    for (ki, &col) in order.iter().take(k).enumerate() {
        for j in 0..n {
            tapers[[ki, j]] = vecs[[j, col]];
        }
    }

    // 4. scipy-compatible sign convention.
    apply_sign_convention(&mut tapers);

    // 5. Concentration ratios via direct quadrature against the sinc kernel.
    let ratios = concentration_ratios(&tapers, nw);

    Ok((tapers, ratios))
}

/// Symmetric tridiagonal eigendecomposition (LAPACK DSTEQR-style).
///
/// `diag` (length N) and `off` (length N-1) define a symmetric tridiagonal
/// matrix. Returns `(eigenvalues, eigenvectors)` where eigenvectors are
/// columns of the returned matrix and eigenvalues are sorted ascending.
fn symmetric_tridiagonal_eigh(
    diag: &[f64],
    off: &[f64],
) -> Result<(Array1<f64>, Array2<f64>), DpssError> {
    let n = diag.len();
    assert_eq!(off.len(), n - 1);

    // Convention: e[i] is the off-diagonal between row i and row i+1, for
    // i ∈ 0..n-1. We allocate length n with e[n-1] = 0 as a sentinel so
    // searching for "negligible e" naturally terminates at the matrix edge.
    let mut d = diag.to_vec();
    let mut e = vec![0.0_f64; n];
    for i in 0..n - 1 {
        e[i] = off[i];
    }

    let mut z = Array2::<f64>::eye(n);
    const MAX_ITERS: usize = 60;

    for l in 0..n {
        let mut iter = 0_usize;
        loop {
            // Find smallest m >= l such that e[m] is negligible relative to
            // its neighbouring diagonal entries. m == n-1 is allowed: e[n-1]
            // is the 0 sentinel, so the trailing eigenvalue deflates.
            let mut m = l;
            while m < n - 1 {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break; // eigenvalue at d[l] converged
            }
            if iter == MAX_ITERS {
                return Err(DpssError::NoConvergence(MAX_ITERS));
            }
            iter += 1;

            // Wilkinson shift from the trailing 2×2 block [d[l], e[l]; e[l], d[l+1]].
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + r.copysign(g));

            let mut s = 1.0_f64;
            let mut c = 1.0_f64;
            let mut p = 0.0_f64;

            // QL sweep from m-1 down to l, chasing the bulge. The Givens
            // rotation at step i operates on rows/columns i and i+1.
            let mut early_break = false;
            let mut i = m;
            while i > l {
                i -= 1;
                let f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    early_break = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                // Accumulate Givens rotation into Z.
                for kk in 0..n {
                    let zf = z[[kk, i + 1]];
                    z[[kk, i + 1]] = s * z[[kk, i]] + c * zf;
                    z[[kk, i]] = c * z[[kk, i]] - s * zf;
                }
            }
            if !early_break {
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        }
    }

    // Sort eigenpairs ascending. Selection sort is fine — N is small.
    for i in 0..n {
        let mut jmin = i;
        for j in (i + 1)..n {
            if d[j] < d[jmin] {
                jmin = j;
            }
        }
        if jmin != i {
            d.swap(i, jmin);
            for row in 0..n {
                let tmp = z[[row, i]];
                z[[row, i]] = z[[row, jmin]];
                z[[row, jmin]] = tmp;
            }
        }
    }

    Ok((Array1::from(d), z))
}

/// scipy-compatible sign normalization. Rows of `tapers` (shape (K, N)) are
/// negated in place so that:
///   * Even k: sum(taper) ≥ 0
///   * Odd  k: the first sample with |w| > max(1e-7, 1/N) is positive
fn apply_sign_convention(tapers: &mut Array2<f64>) {
    let k = tapers.nrows();
    let n = tapers.ncols();
    let thresh = (1.0_f64 / n as f64).max(1e-7);

    for ki in 0..k {
        if ki % 2 == 0 {
            let s: f64 = (0..n).map(|j| tapers[[ki, j]]).sum();
            if s < 0.0 {
                for j in 0..n {
                    tapers[[ki, j]] = -tapers[[ki, j]];
                }
            }
        } else {
            let mut sign_value = 0.0;
            for j in 0..n {
                let v = tapers[[ki, j]];
                if v.abs() > thresh {
                    sign_value = v;
                    break;
                }
            }
            if sign_value < 0.0 {
                for j in 0..n {
                    tapers[[ki, j]] = -tapers[[ki, j]];
                }
            }
        }
    }
}

/// Concentration ratios. For each taper `v_k`,
///   ratio_k = sum_{i,j} v_k[i] · K(i - j) · v_k[j]
/// where K is the prolate kernel:
///   K(m) = 2W · sinc(2W · m)          (m ≠ 0)
///   K(0) = 2W
/// with W = NW/N and `sinc(x) = sin(πx)/(πx)`. Equivalent to the standard
/// Slepian eigenvalue λ_k = ∫_{-W}^{W} |V_k(f)|² df. Direct O(K·N²)
/// quadrature; for N ≤ 1000 this is microseconds.
fn concentration_ratios(tapers: &Array2<f64>, nw: f64) -> Array1<f64> {
    let k = tapers.nrows();
    let n = tapers.ncols();
    let w = nw / n as f64;
    let two_w = 2.0 * w;

    // Prolate kernel values for offset m = 0..n-1.
    let mut kernel = vec![0.0_f64; n];
    kernel[0] = two_w;
    for m in 1..n {
        let arg = std::f64::consts::PI * two_w * m as f64;
        kernel[m] = two_w * (arg.sin() / arg);
    }

    let mut ratios = Array1::<f64>::zeros(k);
    for ki in 0..k {
        let mut acc = 0.0_f64;
        for i in 0..n {
            let vi = tapers[[ki, i]];
            for j in 0..n {
                let m = if i >= j { i - j } else { j - i };
                acc += vi * kernel[m] * tapers[[ki, j]];
            }
        }
        ratios[ki] = acc;
    }
    ratios
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_orthonormal(tapers: &Array2<f64>, tol: f64) {
        let k = tapers.nrows();
        for a in 0..k {
            for b in 0..k {
                let mut dot = 0.0;
                for j in 0..tapers.ncols() {
                    dot += tapers[[a, j]] * tapers[[b, j]];
                }
                let expected = if a == b { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < tol,
                    "tapers {} and {} not orthonormal: dot={}, expected={}",
                    a, b, dot, expected
                );
            }
        }
    }

    fn tridiag_residual(diag: &[f64], off: &[f64], v: &[f64]) -> (f64, f64) {
        let n = diag.len();
        let mut tv = vec![0.0; n];
        for i in 0..n {
            tv[i] = diag[i] * v[i];
            if i > 0 {
                tv[i] += off[i - 1] * v[i - 1];
            }
            if i < n - 1 {
                tv[i] += off[i] * v[i + 1];
            }
        }
        let lambda: f64 = (0..n).map(|i| v[i] * tv[i]).sum();
        let resid: f64 = (0..n)
            .map(|i| (tv[i] - lambda * v[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        (lambda, resid)
    }

    #[test]
    fn rejects_bad_params() {
        assert!(dpss(1, 2.0, 1).is_err());
        assert!(dpss(10, 0.0, 1).is_err());
        assert!(dpss(10, 5.0, 1).is_err()); // NW must be < N/2
        assert!(dpss(10, 2.0, 0).is_err());
        assert!(dpss(10, 2.0, 11).is_err());
    }

    #[test]
    fn small_dpss_is_orthonormal() {
        let (tapers, ratios) = dpss(20, 2.0, 4).unwrap();
        assert_eq!(tapers.shape(), &[4, 20]);
        assert_eq!(ratios.len(), 4);
        assert_orthonormal(&tapers, 1e-10);
    }

    #[test]
    fn dpss_eigenvectors_satisfy_tridiag_equation() {
        let n = 50;
        let nw = 2.5;
        let k = 4;
        let (tapers, _ratios) = dpss(n, nw, k).unwrap();

        let nf = n as f64;
        let cw = (2.0 * std::f64::consts::PI * nw / nf).cos();
        let diag: Vec<f64> = (0..n)
            .map(|i| {
                let m = (nf - 1.0 - 2.0 * i as f64) / 2.0;
                m * m * cw
            })
            .collect();
        let off: Vec<f64> = (0..n - 1)
            .map(|i| {
                let ip = i as f64 + 1.0;
                ip * (nf - ip) / 2.0
            })
            .collect();

        for ki in 0..k {
            let v: Vec<f64> = (0..n).map(|j| tapers[[ki, j]]).collect();
            let (_lambda, resid) = tridiag_residual(&diag, &off, &v);
            assert!(resid < 1e-8, "taper {} residual {} too large", ki, resid);
        }
    }

    #[test]
    fn concentration_ratios_are_descending_in_unit_interval() {
        let (_, ratios) = dpss(100, 3.0, 5).unwrap();
        for ki in 0..ratios.len() {
            assert!(
                ratios[ki] > 0.0 && ratios[ki] <= 1.0 + 1e-9,
                "ratio[{}] = {} out of (0, 1]", ki, ratios[ki]
            );
        }
        for ki in 1..ratios.len() {
            assert!(
                ratios[ki] <= ratios[ki - 1] + 1e-9,
                "ratios not monotone: {} > {} at index {}",
                ratios[ki], ratios[ki - 1], ki
            );
        }
        assert!(ratios[0] > 0.999, "ratio[0] = {} should be ~1.0", ratios[0]);
        assert!(ratios[3] > 0.95, "ratio[3] = {} should still be high", ratios[3]);
    }

    #[test]
    fn sign_convention_even_taper_sum_nonneg() {
        let (tapers, _) = dpss(50, 2.0, 4).unwrap();
        for ki in (0..4).step_by(2) {
            let s: f64 = (0..50).map(|j| tapers[[ki, j]]).sum();
            assert!(s >= -1e-12, "even taper {} has negative sum {}", ki, s);
        }
    }

    #[test]
    fn sign_convention_odd_taper_first_nonsmall_positive() {
        let (tapers, _) = dpss(50, 2.0, 4).unwrap();
        let thresh = (1.0_f64 / 50.0).max(1e-7);
        for ki in (1..4).step_by(2) {
            let mut first = 0.0;
            for j in 0..50 {
                if tapers[[ki, j]].abs() > thresh {
                    first = tapers[[ki, j]];
                    break;
                }
            }
            assert!(
                first > 0.0,
                "odd taper {}: first significant sample = {}",
                ki, first
            );
        }
    }

    #[test]
    fn dynamo_pass1_window_works() {
        // DYNAM-O pass-1 default: 1 s window at 100 Hz → N=100, NW=2, K=3.
        // For NW=2 the "well-concentrated" rule of thumb is 2NW-1 = 3 tapers;
        // the first two sit at ratio ~0.9999, the third typically 0.95–0.99.
        let (tapers, ratios) = dpss(100, 2.0, 3).unwrap();
        assert_eq!(tapers.shape(), &[3, 100]);
        assert_orthonormal(&tapers, 1e-10);
        assert!(ratios[0] > 0.999, "ratio[0] = {}", ratios[0]);
        assert!(ratios[1] > 0.99, "ratio[1] = {}", ratios[1]);
        assert!(ratios[2] > 0.95, "ratio[2] = {}", ratios[2]);
    }
}
