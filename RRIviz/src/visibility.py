# src/visibility.py

import numpy as np
from astropy.coordinates import AltAz, SkyCoord
import astropy.units as u

# Gaussian beamwidth parameter
theta_HPBW = np.deg2rad(10.9)

### Optimized Visibility Calculation Function ###
def calculate_visibility_optimized(antennas, baselines, sources, location, obstime, wavelengths, freqs, sp_index):
    """
    Optimized calculation of complex visibility for each baseline at different frequencies.

    Parameters:
    antennas (dict): Dictionary of antenna positions
    baselines (dict): Dictionary of baselines between antennas
    sources (list): List of source coordinates and fluxes
    location (EarthLocation): Observer's location
    obstime (Time): Observation time
    wavelengths (Quantity): Wavelength array
    freqs (ndarray): Frequency array
    sp_index (float): Spectral index for flux adjustment

    Returns:
    dict: Dictionary of visibilities for each baseline
    """
    visibilities = {key: np.zeros(len(wavelengths), dtype=complex) for key in baselines.keys()}

    # Convert source list to arrays for vectorized operations
    source_coords = SkyCoord([s['coords'] for s in sources])
    source_fluxes = np.array([s['flux'] for s in sources])

    # Transform source coordinates to AltAz frame
    altaz = source_coords.transform_to(AltAz(obstime=obstime, location=location))
    az = altaz.az.rad
    alt = altaz.alt.rad

    # Filter out sources below the horizon
    above_horizon = alt > 0
    az = az[above_horizon]
    alt = alt[above_horizon]
    source_fluxes = source_fluxes[above_horizon]

    theta = np.pi / 2 - alt
    l = np.cos(alt) * np.cos(az)
    m = np.cos(alt) * np.sin(az)
    n = np.sin(alt)

    # Gaussian primary beam pattern
    A_theta = np.exp(- (theta / (np.sqrt(2) * theta_HPBW))**2 ) ** 2

    for i, (wavelength, freq) in enumerate(zip(wavelengths, freqs)):
        # Flux adjusted by spectral index
        flux_adj = source_fluxes * (freq / 50e6)**sp_index

        for (ant1, ant2), baseline in baselines.items():
            u, v, w = np.array(baseline) / wavelength.value  # Convert baseline to units of wavelength

            # Projected baseline component in the direction of the sources
            b_dot_s = u * l + v * m + w * n  # Array over sources

            # Compute the phase term
            phase = np.exp(-2j * np.pi * b_dot_s)

            # Sum over all sources to get the visibility
            visibility = np.sum(flux_adj * A_theta * phase)

            visibilities[(ant1, ant2)][i] = visibility

    return visibilities


### Original Visibility Calculation Function ###
def calculate_visibility_original(antennas, baselines, sources, location, obstime, wavelengths, freqs, sp_index):
    """
    Original calculation of complex visibility for each baseline at different frequencies.

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

### Polarized Visibility Calculation Function ###
def calculate_polarized_visibility(antennas, baselines, sources, location, obstime, wavelengths, freqs, sp_index):
    """
    Polarized calculation of complex visibility for each baseline at different frequencies.
    Takes into account Ex and Ey components.

    Parameters:
    antennas (dict): Dictionary of antenna positions
    baselines (dict): Dictionary of baselines between antennas
    sources (list): List of source coordinates, fluxes, and polarization states
    location (EarthLocation): Observer's location
    obstime (Time): Observation time
    wavelengths (Quantity): Wavelength array
    freqs (ndarray): Frequency array
    sp_index (float): Spectral index for flux adjustment

    Returns:
    dict: Dictionary of visibilities for each baseline with polarization.
    """
    visibilities = {key: np.zeros((len(wavelengths), 2), dtype=complex) for key in baselines.keys()}  # Now storing (Ex, Ey)

    # Convert source list to arrays for vectorized operations
    source_coords = SkyCoord([s['coords'] for s in sources])
    source_fluxes = np.array([s['flux'] for s in sources])
    source_polarization = np.array([s['polarization'] for s in sources])  # Assume polarization (Ex, Ey) provided in source data

    # Transform source coordinates to AltAz frame
    altaz = source_coords.transform_to(AltAz(obstime=obstime, location=location))
    az = altaz.az.rad
    alt = altaz.alt.rad

    # Filter out sources below the horizon
    above_horizon = alt > 0
    az = az[above_horizon]
    alt = alt[above_horizon]
    source_fluxes = source_fluxes[above_horizon]
    source_polarization = source_polarization[above_horizon]

    theta = np.pi / 2 - alt
    l = np.cos(alt) * np.cos(az)
    m = np.cos(alt) * np.sin(az)
    n = np.sin(alt)

    # Gaussian primary beam pattern
    A_theta = np.exp(- (theta / (np.sqrt(2) * theta_HPBW))**2 ) ** 2

    for i, (wavelength, freq) in enumerate(zip(wavelengths, freqs)):
        # Flux adjusted by spectral index
        flux_adj = source_fluxes * (freq / 50e6)**sp_index

        for (ant1, ant2), baseline in baselines.items():
            u, v, w = np.array(baseline) / wavelength.value  # Convert baseline to units of wavelength

            # Projected baseline component in the direction of the sources
            b_dot_s = u * l + v * m + w * n  # Array over sources

            # Compute the phase term
            phase = np.exp(-2j * np.pi * b_dot_s)

            # Calculate Ex and Ey contributions to visibility
            Ex_contrib = np.sum(flux_adj * source_polarization[:, 0] * A_theta * phase)  # Polarization X-component
            Ey_contrib = np.sum(flux_adj * source_polarization[:, 1] * A_theta * phase)  # Polarization Y-component

            # Store both Ex and Ey visibilities
            visibilities[(ant1, ant2)][i, 0] = Ex_contrib
            visibilities[(ant1, ant2)][i, 1] = Ey_contrib

    return visibilities


### Common Function to Calculate Modulus and Phase ###
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
