# src/skymodel.py

import pygdsm
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
import astropy.units as u


def collect_sky_model_for_time(location, obstime, frequency=76, fov_radius_deg=5, gleam_sources=None):
    """
    Generates the Global Sky Model data at a given frequency and time,
    and returns a function to execute the plot later.

    Parameters:
    - location (EarthLocation): Observer's location.
    - obstime (Time): Observation time.
    - frequency (float): Frequency in MHz for the sky model.
    - fov_radius_deg (float): Radius of the field of view in degrees.
    - gleam_sources (list): List of GLEAM sources to overlay on the plot.

    Returns:
    - A function that, when called, will plot the sky model for the specified time.
    """
    # Generate the Global Sky Model at the specified frequency
    gsm_2008 = pygdsm.GlobalSkyModel(freq_unit='MHz')
    sky_map = gsm_2008.generate(frequency)

    # Determine the number of pixels (nside parameter) for the HEALPix map
    nside = 256  # Higher resolution

    # Generate the HEALPix map
    healpix_map = hp.ud_grade(sky_map, nside_out=nside)

    # Rotate the HEALPix map to Equatorial coordinates
    rotator = hp.Rotator(coord=['G', 'C'])
    equatorial_map = rotator.rotate_map_pixel(healpix_map)

    # HERA declination
    hera_dec = location.lat.deg
    fov_radius_deg = fov_radius_deg  # Field of view radius

    def plot_sky_model():
        # Create a new figure for each plot
        fig = plt.figure(figsize=(14, 7))

        # Plot the equatorial sky model in Cartesian coordinates
        hp.cartview(
            equatorial_map,
            coord=['C'],
            title=f"Time: {obstime.iso} | GSM Sky Model at {frequency} MHz",
            norm='log',
            xsize=2000,
            flip='astro',
            cmap='inferno',
            return_projected_map=False
        )

        # Define the zenith in altaz coordinates at the given observation time
        zenith = SkyCoord(
            alt=90 * u.deg,
            az=0 * u.deg,
            frame=AltAz(obstime=obstime, location=location)
        )

        # Convert the zenith to equatorial coordinates
        zenith_radec = zenith.transform_to('icrs')
        ra_highlight = zenith_radec.ra.deg  # Right Ascension of the zenith

        # Calculate the RA and Dec boundaries for the highlighted area
        ra_range_highlight = np.linspace(
            ra_highlight - fov_radius_deg, ra_highlight + fov_radius_deg, 100
        )
        dec_upper_highlight = [hera_dec + fov_radius_deg] * 100
        dec_lower_highlight = [hera_dec - fov_radius_deg] * 100

        # Plot the upper, lower, and side boundaries of the highlighted area
        hp.visufunc.projplot(
            ra_range_highlight,
            dec_upper_highlight,
            lonlat=True,
            color='red',
            linewidth=1.5,
            alpha=1.0
        )
        hp.visufunc.projplot(
            ra_range_highlight,
            dec_lower_highlight,
            lonlat=True,
            color='red',
            linewidth=1.5,
            alpha=1.0
        )
        hp.visufunc.projplot(
            [ra_highlight - fov_radius_deg, ra_highlight - fov_radius_deg],
            [hera_dec - fov_radius_deg, hera_dec + fov_radius_deg],
            lonlat=True,
            color='red',
            linewidth=1.5,
            alpha=1.0
        )
        hp.visufunc.projplot(
            [ra_highlight + fov_radius_deg, ra_highlight + fov_radius_deg],
            [hera_dec - fov_radius_deg, hera_dec + fov_radius_deg],
            lonlat=True,
            color='red',
            linewidth=1.5,
            alpha=1.0
        )

        # If GLEAM sources are provided, overlay them on the map and label them
        if gleam_sources:
            ra_gleam = [src['coords'].ra.deg for src in gleam_sources]
            dec_gleam = [src['coords'].dec.deg for src in gleam_sources]
            flux_gleam = [src['flux'] for src in gleam_sources]

            # Plot GLEAM sources as scatter points, size scaled by flux
            hp.visufunc.projscatter(
                ra_gleam,
                dec_gleam,
                lonlat=True,
                s=np.log10(flux_gleam) * 10,  # Scale size by log of flux
                color='blue',
                alpha=0.7,
                label='GLEAM sources'
            )

            # Label the GLEAM sources with flux
            for ra, dec, flux in zip(ra_gleam, dec_gleam, flux_gleam):
                plt.annotate( 
                    f"{flux:.1f} Jy",  # Label with flux
                    xy=(ra, dec),
                    xycoords='data',
                    xytext=(5, 5),  # Offset label slightly
                    textcoords='offset points',
                    fontsize=8,
                    color='white'
                )

        # Return the figure object to be used later
        return fig

    return plot_sky_model
