# tests/test_observation.py

import unittest
from src.observation import get_location_and_time
import astropy.units as u

class TestObservation(unittest.TestCase):

    def test_default_location(self):
        location, obstime = get_location_and_time()
        self.assertAlmostEqual(location.lat.deg, -30.72152777777791)
        self.assertAlmostEqual(location.lon.deg, 21.428305555555557)
        self.assertAlmostEqual(location.height.to(u.m).value, 1073.0)

    def test_custom_location(self):
        lat = -30.0
        lon = 21.0
        height = 1000.0
        location, obstime = get_location_and_time(lat, lon, height)
        self.assertAlmostEqual(location.lat.deg, lat)
        self.assertAlmostEqual(location.lon.deg, lon)
        self.assertAlmostEqual(location.height.to(u.m).value, height)

if __name__ == '__main__':
    unittest.main()
