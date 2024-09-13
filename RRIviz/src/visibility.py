# src/visibility.py

import numpy as np
from astropy.coordinates import AltAz

# Gaussian beamwidth parameter
theta_HPBW = np.deg2rad(10.9)

def calculate_visibility(antennas, baselines, sources, location, obstime, wavelengths, freqs, sp_index):
    """
    Calculate the complex visibility for each baseline at different frequencies.

    Parameters:
    antennas (dict): Dictionary of antenna positions
    baselines (dict): Dictionary of baselines between antennas
    sources (list): List of source coordinates and fluxes
    location (EarthLocation): Observer's location
    obstime (Time): Observation time
    wavelengths (Quantity): Wavelength array
    freqs (ndarray): Frequency array
    sp_index (int): Spectral index for flux adjustment

    Returns:
    dict: Dictionary of visibilities for each baseline
    """
    visibilities = {key: np.zeros(len(wavelengths), dtype=complex) for key in baselines.keys()}
    
    for i, (wavelength, freq) in enumerate(zip(wavelengths, freqs)):
        for (ant1, ant2), baseline in baselines.items():
            u, v, w = np.array(baseline) / wavelength.value  # Convert baseline to units of wavelength
            
            visibility = 0j
            for source in sources:
                altaz = source['coords'].transform_to(AltAz(obstime=obstime, location=location))
                az, alt = altaz.az.rad, altaz.alt.rad
                
                theta = np.deg2rad(90 - altaz.alt.deg)
                l = np.cos(alt) * np.cos(az)
                m = np.cos(alt) * np.sin(az)
                n = np.sin(alt)
                
                # Projected baseline component in the direction of the source
                b_dot_s = u * l + v * m + w * n
                
                # Gaussian primary beam pattern
                A_theta = (np.exp(- (theta / (np.sqrt(2) * theta_HPBW))**2))**2
                
                # Flux adjusted by spectral index
                flux_adj = source['flux'] * (freq / 50e6)**sp_index
                
                # Add contribution to visibility with beam pattern and flux adjustment
                visibility += flux_adj * A_theta * np.exp(-2j * np.pi * b_dot_s)
            
            visibilities[(ant1, ant2)][i] = visibility
    
    return visibilities

def calculate_modulus_phase(visibilities):
    """
    Calculate the modulus and phase of visibilities.

    Parameters:
    visibilities (dict): Dictionary of visibilities for each baseline

    Returns:
    tuple: Two dictionaries, one for modulus and one for phase of visibilities
    """
    moduli = {key: np.abs(val) for key, val in visibilities.items()}
    phases = {key: np.angle(val) for key, val in visibilities.items()}
    return moduli, phases
