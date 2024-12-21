# tests/test_antenna.py

import unittest
from unittest.mock import mock_open, patch
from src.antenna import read_antenna_positions

class TestAntenna(unittest.TestCase):

    def setUp(self):
        # Correct sample antenna position file content
        self.correct_file_content = """Name  Number   BeamID   E          N          U         
HH136      136        0  -156.5976     2.9439    -0.1819
HH140      140        0   -98.1662     3.1671    -0.3008
HH121      121        0   -90.8139    -9.4618    -0.1707
"""

        # Malformed file content (missing a value)
        self.malformed_file_content = """Name  Number   BeamID   E          N          U         
HH136      136        0  -156.5976     2.9439
"""

        # Empty file content
        self.empty_file_content = ""

        # File content with extra invalid lines
        self.extra_invalid_content = """Name  Number   BeamID   E          N          U         
HH136      136        0  -156.5976     2.9439    -0.1819
INVALID LINE HERE
HH140      140        0   -98.1662     3.1671    -0.3008
"""

    def test_read_antenna_positions_correct(self):
        with patch("builtins.open", mock_open(read_data=self.correct_file_content)) as mock_file:
            antennas = read_antenna_positions("dummy_path.txt")
            expected = {
                0: {"Name": "HH136", "Number": 136, "BeamID": 0, "Position": (-156.5976, 2.9439, -0.1819)},
                1: {"Name": "HH140", "Number": 140, "BeamID": 0, "Position": (-98.1662, 3.1671, -0.3008)},
                2: {"Name": "HH121", "Number": 121, "BeamID": 0, "Position": (-90.8139, -9.4618, -0.1707)},
            }

            for key, value in expected.items():
                with self.subTest(antenna_index=key):
                    self.assertEqual(
                        antennas[key],
                        value,
                        f"Failed to parse antenna {key}. Expected {value}, got {antennas[key]}"
                    )

            mock_file.assert_called_once_with("dummy_path.txt", "r")

    def test_read_antenna_positions_empty(self):
        with patch("builtins.open", mock_open(read_data=self.empty_file_content)) as mock_file:
            with self.assertRaises(ValueError, msg="Expected ValueError for an empty file but none was raised."):
                read_antenna_positions("dummy_path.txt")

    def test_read_antenna_positions_malformed(self):
        with patch("builtins.open", mock_open(read_data=self.malformed_file_content)) as mock_file:
            with self.assertRaises(ValueError, msg="Expected ValueError for a malformed file but none was raised."):
                read_antenna_positions("dummy_path.txt")

    def test_read_antenna_positions_extra_invalid(self):
        with patch("builtins.open", mock_open(read_data=self.extra_invalid_content)) as mock_file:
            with self.assertRaises(ValueError, msg="Expected ValueError for a file with extra invalid lines but none was raised."):
                read_antenna_positions("dummy_path.txt")


if __name__ == "__main__":
    unittest.main()

