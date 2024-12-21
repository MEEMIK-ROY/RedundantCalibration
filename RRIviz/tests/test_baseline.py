import unittest
import numpy as np
from src.baseline import generate_baselines

class TestBaseline(unittest.TestCase):

    def setUp(self):
        # Sample valid antenna data
        self.antennas = {
            0: {"Name": "HH136", "Number": 136, "BeamID": 0, "Position": (-156.5976, 2.9439, -0.1819)},
            1: {"Name": "HH140", "Number": 140, "BeamID": 0, "Position": (-98.1662, 3.1671, -0.3008)},
            2: {"Name": "HH121", "Number": 121, "BeamID": 0, "Position": (-90.8139, -9.4618, -0.1707)},
        }

        # Malformed antenna data (missing Number or Position)
        self.malformed_antennas = {
            0: {"Name": "HH136", "BeamID": 0},  # Missing Number and Position
            1: {"Name": "HH140", "Number": 140, "BeamID": 0},  # Missing Position
        }

        # Antennas with empty positions
        self.empty_antennas = {}

    def test_generate_baselines_valid(self):
        """Test baseline generation with valid antenna data."""
        baselines = generate_baselines(self.antennas)
        expected = {
            (136, 136): np.array([0.0, 0.0, 0.0]),
            (136, 140): np.array([58.4314, 0.2232, -0.1189]),
            (121, 136): np.array([-65.7837, 12.4057, -0.0112]),
            (140, 140): np.array([0.0, 0.0, 0.0]),
            (121, 140): np.array([-7.3523, 12.6289, -0.1301]),
            (121, 121): np.array([0.0, 0.0, 0.0]),
        }

        # Verify each expected baseline
        for key, value in expected.items():
            with self.subTest(baseline=key):
                np.testing.assert_array_almost_equal(
                    baselines[key],
                    value,
                    decimal=4,
                    err_msg=f"Baseline {key} mismatch. Expected {value}, got {baselines[key]}."
                )

        # Verify that no additional baselines were generated
        self.assertEqual(
            set(baselines.keys()),
            set(expected.keys()),
            f"Unexpected baselines generated. Expected keys: {set(expected.keys())}, but got: {set(baselines.keys())}"
        )

    def test_generate_baselines_empty(self):
        """Test baseline generation with no antennas."""
        baselines = generate_baselines(self.empty_antennas)
        self.assertEqual(
            baselines,
            {},
            "Expected no baselines for empty antenna data but got some baselines."
        )

    def test_generate_baselines_malformed(self):
        """Test baseline generation with malformed antenna data."""
        with self.assertRaises(KeyError, msg="Expected KeyError for missing 'Number' or 'Position' fields in antenna data."):
            generate_baselines(self.malformed_antennas)


if __name__ == "__main__":
    unittest.main()
