"""Regression tests for NaN-containing pure-Python multitaper windows."""

import os
import sys
import unittest

import numpy as np
from scipy.signal.windows import dpss

# Run directly from the repository without installing a package.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from multitaper_spectrogram_python import (
    calc_mts_batch,
    calc_mts_segment,
    multitaper_spectrogram,
)


class BatchedNanWindowTests(unittest.TestCase):
    def test_batch_matches_single_window_nan_contract(self):
        rng = np.random.default_rng(42)
        data_segments = rng.standard_normal((8, 16))
        data_segments[1, 4] = np.nan
        data_segments[2:4, :] = np.nan
        data_segments[4, 7] = np.nan
        data_segments[6, :] = 0

        num_tapers = 3
        nfft = 32
        tapers, eigen = dpss(
            16, 2, num_tapers, return_ratios=True
        )
        tapers_t = tapers.T
        eigen = eigen.reshape(num_tapers, 1)
        freq_inds = np.arange(nfft // 2 + 1) % 3 != 2

        for detrend_opt in ("linear", "constant", "off"):
            for weighting in ("unity", "eigen"):
                with self.subTest(
                    detrend_opt=detrend_opt, weighting=weighting
                ):
                    if weighting == "unity":
                        weights = np.ones((num_tapers, 1)) / num_tapers
                    else:
                        weights = eigen / num_tapers

                    params = (
                        tapers_t,
                        nfft,
                        freq_inds,
                        detrend_opt,
                        num_tapers,
                        eigen,
                        weighting,
                        weights,
                    )
                    expected = np.vstack([
                        calc_mts_segment(segment, *params)
                        for segment in data_segments
                    ])
                    actual = calc_mts_batch(
                        data_segments, *params, batch_size=2
                    )

                    np.testing.assert_allclose(
                        actual,
                        expected,
                        rtol=1e-12,
                        atol=1e-12,
                        equal_nan=True,
                    )
                    self.assertTrue(np.isnan(actual[1:5]).all())
                    self.assertTrue((actual[6] == 0).all())

    def test_full_python_backend_invalidates_overlapping_nan_windows(self):
        fs = 16.0
        time = np.arange(64) / fs
        data = np.sin(2 * np.pi * 2 * time)
        data[20] = np.nan
        data[48:] = 0

        spect, _, _ = multitaper_spectrogram(
            data,
            fs,
            frequency_range=[0, fs / 2],
            time_bandwidth=2,
            num_tapers=3,
            window_params=[1, 0.5],
            min_nfft=0,
            detrend_opt="linear",
            weighting="unity",
            plot_on=False,
            verbose=False,
            use_rust=False,
        )

        expected_nan_columns = np.array([
            False, True, True, False, False, False, False
        ])
        np.testing.assert_array_equal(
            np.isnan(spect).all(axis=0), expected_nan_columns
        )
        self.assertTrue(np.isfinite(spect[:, ~expected_nan_columns]).all())
        self.assertTrue((spect[:, -1] == 0).all())


if __name__ == "__main__":
    unittest.main()
