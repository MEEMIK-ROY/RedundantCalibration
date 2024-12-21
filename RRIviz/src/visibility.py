# src/visibility.py

import numpy as np
import sys
from astropy.coordinates import AltAz, SkyCoord
import astropy.units as u


def calculate_A_theta(theta, theta_HPBW):
    """
    Calculate the Gaussian primary beam pattern A(theta).

    Parameters:
    theta (ndarray): Zenith angles in radians.
    theta_HPBW (float): Half Power Beam Width (HPBW) in radians.

    Returns:
    ndarray: The primary beam pattern A(theta) evaluated at theta.
    """
    A_theta = np.exp(-((theta / (np.sqrt(2) * theta_HPBW)) ** 2)) ** 2
    return A_theta


### Optimized Visibility Calculation Function ###


def calculate_visibility_optimized(
    antennas, baselines, sources, location, obstime, wavelengths, freqs, theta_HPBW
):
    """
    Optimized calculation of complex visibility for each baseline at different frequencies.

    Parameters:
    antennas (dict): Dictionary of antenna positions.
    baselines (dict): Dictionary of baselines between antennas.
    sources (list): List of source dictionaries containing 'coords', 'flux', and 'spectral_index'.
    location (EarthLocation): Observer's geographical location.
    obstime (Time): Observation time.
    wavelengths (Quantity): Wavelength array corresponding to the frequencies.
    freqs (ndarray): Frequency array in Hz.
    theta_HPBW (float): Half Power Beam Width (HPBW) in radians.

    Returns:
    dict: Dictionary of visibilities for each baseline as arrays over frequencies.
    """
    # Initialize visibilities dictionary with zeros for each baseline and frequency
    visibilities = {
        key: np.zeros(len(wavelengths), dtype=complex) for key in baselines.keys()
    }

    # Handle the case where there are no sources
    if not sources:
        return visibilities

    # Convert source list to arrays for vectorized operations
    source_coords = SkyCoord([s["coords"] for s in sources])
    source_fluxes = np.array([s["flux"] for s in sources])
    source_spectral_indices = np.array([s["spectral_index"] for s in sources])

    # Reference frequency for spectral index calculation is 76 MHz
    reference_freq = 76e6  # in Hz

    # Transform source coordinates to AltAz frame at the observation time and location
    altaz = source_coords.transform_to(AltAz(obstime=obstime, location=location))
    az = altaz.az.rad  # Azimuth in radians
    alt = altaz.alt.rad  # Altitude in radians

    # Filter out sources below the horizon (altitude <= 0)
    above_horizon = alt > 0
    az = az[above_horizon]
    alt = alt[above_horizon]
    source_fluxes = source_fluxes[above_horizon]
    source_spectral_indices = source_spectral_indices[above_horizon]

    # Calculate zenith angle theta
    theta = np.pi / 2 - alt
    # Direction cosines for the sources
    l = np.cos(alt) * np.sin(az)
    m = np.cos(alt) * np.cos(az)
    n = np.sin(alt)

    # Gaussian primary beam pattern as a function of theta
    A_theta = calculate_A_theta(theta, theta_HPBW)

    # Loop over each frequency and wavelength
    for i, (wavelength, freq) in enumerate(zip(wavelengths, freqs)):
        # Adjust fluxes by the spectral indices at the current frequency
        flux_adj = source_fluxes * (freq / reference_freq) ** source_spectral_indices

        # Loop over each baseline
        for (ant1, ant2), baseline in baselines.items():
            # Convert baseline vector to units of wavelength
            u, v, w = np.array(baseline) / wavelength.value

            # Projected baseline component in the direction of the sources
            b_dot_s = u * l + v * m + w * n  # Array over sources

            # Compute the phase term e^{-2πi * (b ⋅ s)}
            phase = np.exp(-2j * np.pi * b_dot_s)

            # Sum over all sources to get the visibility for this baseline and frequency
            visibility = np.sum(flux_adj * A_theta * phase)

            # Store the visibility
            visibilities[(ant1, ant2)][i] = visibility
    # # Calculate total memory usage
    # total_memory_bytes = sys.getsizeof(visibilities) + sum(
    #     sys.getsizeof(key) + sys.getsizeof(value) + value.nbytes for key, value in visibilities.items()
    # )
    # total_memory_mb = total_memory_bytes / (1024 * 1024)
    # print(f"Total memory used by visibilities: {total_memory_mb:.4f} MB")

    return visibilities


### Original Visibility Calculation Function ###


def calculate_visibility_original(
    antennas, baselines, sources, location, obstime, wavelengths, freqs, theta_HPBW
):
    """
    Original calculation of complex visibility for each baseline at different frequencies.

    Parameters:
    antennas (dict): Dictionary of antenna positions.
    baselines (dict): Dictionary of baselines between antennas.
    sources (list): List of source dictionaries containing 'coords', 'flux', and 'spectral_index'.
    location (EarthLocation): Observer's geographical location.
    obstime (Time): Observation time.
    wavelengths (Quantity): Wavelength array corresponding to the frequencies.
    freqs (ndarray): Frequency array in Hz.
    theta_HPBW (float): Half Power Beam Width (HPBW) in radians.

    Returns:
    dict: Dictionary of visibilities for each baseline as arrays over frequencies.
    """
    # Initialize visibilities dictionary with zeros for each baseline and frequency
    visibilities = {
        key: np.zeros(len(wavelengths), dtype=complex) for key in baselines.keys()
    }

    # Reference frequency for spectral index calculation is 76 MHz
    reference_freq = 76e6  # in Hz

    # Loop over each frequency and wavelength
    for i, (wavelength, freq) in enumerate(zip(wavelengths, freqs)):
        # Loop over each baseline
        for (ant1, ant2), baseline in baselines.items():
            # Convert baseline vector to units of wavelength
            u, v, w = np.array(baseline) / wavelength.value

            visibility = 0j  # Initialize visibility for this baseline and frequency
            # Loop over each source
            for source in sources:
                # Transform source coordinate to AltAz frame at the observation time and location
                altaz = source["coords"].transform_to(
                    AltAz(obstime=obstime, location=location)
                )
                az, alt = altaz.az.rad, altaz.alt.rad  # Azimuth and altitude in radians

                if alt <= 0:
                    continue  # Skip sources below the horizon

                # Calculate zenith angle theta
                theta = np.pi / 2 - alt
                # Direction cosines for the source
                l = np.cos(alt) * np.sin(az)
                m = np.cos(alt) * np.cos(az)
                n = np.sin(alt)

                # Projected baseline component in the direction of the source
                b_dot_s = u * l + v * m + w * n

                # Gaussian primary beam pattern
                A_theta = calculate_A_theta(theta, theta_HPBW)

                # Flux adjusted by the spectral index at the current frequency
                spectral_index = source["spectral_index"]
                flux_adj = source["flux"] * (freq / reference_freq) ** spectral_index

                # Add contribution to visibility with beam pattern and flux adjustment
                visibility += flux_adj * A_theta * np.exp(-2j * np.pi * b_dot_s)

            # Store the visibility
            visibilities[(ant1, ant2)][i] = visibility

    return visibilities


### Polarized Visibility Calculation Function ###


def calculate_polarized_visibility(
    antennas, baselines, sources, location, obstime, wavelengths, freqs, theta_HPBW
):
    """
    Polarized calculation of complex visibility for each baseline at different frequencies.
    Takes into account Ex and Ey components.

    Parameters:
    antennas (dict): Dictionary of antenna positions.
    baselines (dict): Dictionary of baselines between antennas.
    sources (list): List of source dictionaries containing 'coords', 'flux', 'polarization', and 'spectral_index'.
    location (EarthLocation): Observer's geographical location.
    obstime (Time): Observation time.
    wavelengths (Quantity): Wavelength array corresponding to the frequencies.
    freqs (ndarray): Frequency array in Hz.
    theta_HPBW (float): Half Power Beam Width (HPBW) in radians.

    Returns:
    dict: Dictionary of visibilities for each baseline with polarization components.
          Each value is an array of shape (len(wavelengths), 2) for Ex and Ey components.
    """
    # Initialize visibilities dictionary with zeros for each baseline, frequency, and polarization component
    visibilities = {
        key: np.zeros((len(wavelengths), 2), dtype=complex) for key in baselines.keys()
    }  # Now storing (Ex, Ey)

    # Convert source list to arrays for vectorized operations
    source_coords = SkyCoord([s["coords"] for s in sources])
    source_fluxes = np.array([s["flux"] for s in sources])
    source_polarization = np.array(
        [s["polarization"] for s in sources]
    )  # Assume polarization (Ex, Ey) provided in source data
    source_spectral_indices = np.array([s["spectral_index"] for s in sources])

    # Reference frequency for spectral index calculation is 76 MHz
    reference_freq = 76e6  # in Hz

    # Transform source coordinates to AltAz frame at the observation time and location
    altaz = source_coords.transform_to(AltAz(obstime=obstime, location=location))
    az = altaz.az.rad  # Azimuth in radians
    alt = altaz.alt.rad  # Altitude in radians

    # Filter out sources below the horizon (altitude <= 0)
    above_horizon = alt > 0
    az = az[above_horizon]
    alt = alt[above_horizon]
    source_fluxes = source_fluxes[above_horizon]
    source_polarization = source_polarization[above_horizon]
    source_spectral_indices = source_spectral_indices[above_horizon]

    # Calculate zenith angle theta
    theta = np.pi / 2 - alt
    # Direction cosines for the sources
    l = np.cos(alt) * np.cos(az)
    m = np.cos(alt) * np.sin(az)
    n = np.sin(alt)

    # Gaussian primary beam pattern as a function of theta
    A_theta = calculate_A_theta(theta, theta_HPBW)

    # Loop over each frequency and wavelength
    for i, (wavelength, freq) in enumerate(zip(wavelengths, freqs)):
        # Adjust fluxes by individual spectral indices at the current frequency
        flux_adj = source_fluxes * (freq / reference_freq) ** source_spectral_indices

        # Loop over each baseline
        for (ant1, ant2), baseline in baselines.items():
            # Convert baseline vector to units of wavelength
            u, v, w = np.array(baseline) / wavelength.value

            # Projected baseline component in the direction of the sources
            b_dot_s = u * l + v * m + w * n  # Array over sources

            # Compute the phase term e^{-2πi * (b ⋅ s)}
            phase = np.exp(-2j * np.pi * b_dot_s)

            # Calculate Ex and Ey contributions to visibility
            Ex_contrib = np.sum(
                flux_adj * source_polarization[:, 0] * A_theta * phase
            )  # Polarization X-component
            Ey_contrib = np.sum(
                flux_adj * source_polarization[:, 1] * A_theta * phase
            )  # Polarization Y-component

            # Store both Ex and Ey visibilities
            visibilities[(ant1, ant2)][i, 0] = Ex_contrib
            visibilities[(ant1, ant2)][i, 1] = Ey_contrib

    return visibilities


### Common Function to Calculate Modulus and Phase ###


def calculate_modulus_phase(visibilities):
    """
    Calculate the modulus (amplitude) and phase of visibilities.

    Parameters:
    visibilities (dict): Dictionary of visibilities for each baseline.

    Returns:
    tuple: Two dictionaries,
           - moduli: Dictionary of amplitudes of visibilities for each baseline.
           - phases: Dictionary of phases (in radians) of visibilities for each baseline.
    """
    moduli = {key: np.abs(val) for key, val in visibilities.items()}
    phases = {key: np.angle(val) for key, val in visibilities.items()}
    return moduli, phases
