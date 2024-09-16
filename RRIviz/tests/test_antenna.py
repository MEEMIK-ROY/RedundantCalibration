# tests/test_antenna.py

import unittest
from src.antenna import generate_antennas

class TestAntenna(unittest.TestCase):

    def test_generate_default_antennas(self):
        antennas = generate_antennas()
        self.assertEqual(len(antennas), 3)
        self.assertEqual(antennas[0], (0.0, 0, 0))
        self.assertEqual(antennas[1], (14.0, 0, 0))
        self.assertEqual(antennas[2], (28.0, 0, 0))

    def test_generate_custom_spacing(self):
        num_antennas = 4
        spacing = 20.0
        antennas = generate_antennas(num_antennas=num_antennas, spacing=spacing)
        self.assertEqual(len(antennas), num_antennas)
        self.assertEqual(antennas[0], (0.0, 0, 0))
        self.assertEqual(antennas[1], (20.0, 0, 0))
        self.assertEqual(antennas[2], (40.0, 0, 0))
        self.assertEqual(antennas[3], (60.0, 0, 0))

if __name__ == '__main__':
    unittest.main()
