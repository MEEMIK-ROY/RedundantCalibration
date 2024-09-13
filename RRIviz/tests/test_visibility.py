# tests/test_visibility.py

import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from src.visibility import calculate_visibility, calculate_modulus_phase

def test_calculate_visibility():
    # Setting up test data
    antennas = {0: (0, 0, 0), 1: (14, 0, 0), 2: (28, 0, 0)}
    baselines = {(0, 1): np.array([14, 0, 0]), (0, 2): np.array([28, 0, 0]), (1, 2): np.array([14, 0, 0])}
    sources = [{'coords': SkyCoord(ra=0*u.deg, dec=-30.72152777777791*u.deg), 'flux': 2}]
    
    location = EarthLocation(lat=-30.72152777777791*u.deg, lon=21.428305555555557*u.deg, height=1073*u.m)
    obstime = Time("2024-06-05T00:00:00", format="isot", scale="utc")
    wavelengths = np.array([3]) * u.m  # Placeholder for one wavelength
    freqs = np.array([50e6])     # Placeholder for one frequency
    sp_index = 0

    # Run the visibility calculation
    visibilities = calculate_visibility(antennas, baselines, sources, location, obstime, wavelengths, freqs, sp_index)
    
    # Check that visibilities are calculated and have the correct shape
    assert isinstance(visibilities, dict)
    assert len(visibilities) == 3
    for key in visibilities:
        assert len(visibilities[key]) == len(wavelengths)

def test_calculate_modulus_phase():
    # Placeholder visibility data
    visibilities = {(0, 1): np.array([1+1j, 1-1j]), (0, 2): np.array([0+0j]), (1, 2): np.array([2+2j])}

    # Run modulus and phase calculation
    moduli, phases = calculate_modulus_phase(visibilities)

    # Check that modulus and phase values are calculated correctly
    assert np.allclose(moduli[(0, 1)], np.array([np.sqrt(2), np.sqrt(2)]))
    assert np.allclose(phases[(0, 1)], np.array([np.pi/4, -np.pi/4]))
