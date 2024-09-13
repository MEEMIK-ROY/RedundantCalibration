# tests/test_baseline.py

import numpy as np
from src.baseline import baselines

def test_baseline_lengths():
    # Test if the correct number of baselines are present
    assert len(baselines) == 3

def test_baseline_values():
    # Test if the baseline values are correct
    np.testing.assert_array_equal(baselines[(0, 1)], np.array([14, 0, 0]))
    np.testing.assert_array_equal(baselines[(0, 2)], np.array([28, 0, 0]))
    np.testing.assert_array_equal(baselines[(1, 2)], np.array([14, 0, 0]))
