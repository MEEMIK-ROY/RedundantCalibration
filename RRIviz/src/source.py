# src/source.py

from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from pympler import asizeof
import numpy as np
import healpy as hp
from healpy.rotator import Rotator
from pygdsm import GlobalSkyModel
from beams import calculate_A_theta

# Test sources (used when not loading GLEAM catalog)
test_sources = [
    {
        "coords": SkyCoord(ra=0 * u.deg, dec=-30.72152777777791 * u.deg),
        "flux": 2,
        "spectral_index": -0.8,
    },
    {
        "coords": SkyCoord(ra=120 * u.deg, dec=-30.72152777777791 * u.deg),
        "flux": 4,
        "spectral_index": -0.8,
    },
    {
        "coords": SkyCoord(ra=240 * u.deg, dec=-30.72152777777791 * u.deg),
        "flux": 6,
        "spectral_index": -0.8,
    },
]


# GLEAM catalog loading function
def load_gleam_catalog(flux_limit=1.0):
    """
    Loads the GLEAM Extragalactic catalog from VizieR using astroquery.
    Applies a flux limit to reduce the number of sources.

    Parameters:
    flux_limit (float): Minimum flux (in Jy) to include a source.

    Returns:
    list: List of sources with their coordinates, fluxes, and spectral indices.
    """
    # Set up Vizier query
    Vizier.ROW_LIMIT = -1  # No row limit to retrieve all sources

    # Retrieve the GLEAM Extragalactic Catalog
    catalog_list = Vizier.find_catalogs("VIII/100/gleamegc")
    if not catalog_list:
        raise Exception("GLEAM catalog not found in VizieR.")
    catalog = Vizier.get_catalogs("VIII/100/gleamegc")[0]

    # Apply flux limit and extract positions and fluxes
    sources = []
    for row in catalog:
        flux = row["Fp076"]  # Wide-band peak flux density at 76 MHz in Jy
        if flux >= flux_limit:
            ra = row["RAJ2000"] * u.deg
            dec = row["DEJ2000"] * u.deg
            coords = SkyCoord(ra=ra, dec=dec)

            # Retrieve the spectral index from the 'Alpha' column
            if "alpha" in row.colnames and isinstance(row["alpha"], (float, np.float32, np.float64)) and np.isfinite(row["alpha"]):
                spectral_index = row["alpha"]
            else:
                spectral_index = 0.0  # Assign a default value if missing
                
            print(f"Source RA: {ra}, Dec: {dec}, Flux: {flux}, Spectral Index: {spectral_index}")

            # Append source dictionary to the list
            sources.append(
                {"coords": coords, "flux": flux, "spectral_index": spectral_index}
            )

    # Measure the memory size of the loaded GLEAM catalog
    gleam_size = asizeof.asizeof(sources)
    print(
        f"GLEAM catalog loaded with {len(sources)} sources. Total size: {gleam_size / (1024**2):.2f} MB"
    )

    return sources




def load_gsm2008_sources(frequency=76e6, nside=32, flux_limit=1.0):
    """
    Processes the GSM2008 Healpix map and converts it into sources compatible with the app.

    Parameters:
    - frequency (float): Frequency in Hz for flux density calculations.
    - nside (int): Healpix resolution parameter.
    - flux_limit (float): Minimum flux (in Jy) to include a source.

    Returns:
    list: List of sources with "coords", "flux", and "spectral_index" keys.
    """

    gsm = GlobalSkyModel(freq_unit="Hz") # Healpix map in Galactic coordinates with nside=512
    sky_map = gsm.generate(frequency) 

    # Create a rotator to transform from Galactic to Equatorial coordinates
    rot = Rotator(coord=['G', 'C'])
    
    # Downgrade the map
    downgraded_map = hp.ud_grade(sky_map, nside_out=nside)

    # Rotate the map to Equatorial coordinates
    equatorial = rot.rotate_map_pixel(downgraded_map)
    
    
    # Get pixel indices for nside
    npix = hp.nside2npix(nside)
    pixel_indices = np.arange(npix)
    
    # Get theta and phi values for each pixel in galactic coordinates
    theta_gal, phi_gal = hp.pix2ang(nside, pixel_indices)

    # Rotate the pixel to Equatorial coordinates
    theta_eq, phi_eq = rot(theta_gal, phi_gal)
    
    
    # Convert to RA/Dec
    ra = np.rad2deg(phi_eq)
    dec = 90 - np.rad2deg(theta_eq)

    # Calculate flux density
    total_pixels = hp.nside2pixarea(nside)
    solid_angle = (4 * np.pi) / total_pixels  # Solid angle of each pixel in steradians
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    c = 299792458  # Speed of light in m/s
    flux_density = (
        ((2 * k_B * equatorial * (frequency ** 2)) / (c ** 2)) * solid_angle
    )  # Flux density in W/m^2/Hz

    # Convert flux density to Jy (1 Jy = 1e-26 W/m^2/Hz)
    flux_density_jy = flux_density / 1e-26


    # Debugging: Log the range of flux densities
    print(f"Flux density range: {flux_density_jy.min():.2f} Jy to {flux_density_jy.max():.2f} Jy")

    # Apply flux limit condition
    valid_indices = flux_density_jy >= flux_limit
    if valid_indices.sum() == 0:
        print("Warning: No sources meet the flux limit condition.")
    else:
        print(f"{valid_indices.sum()} sources meet the flux limit condition.")

    # Apply flux limit condition
    valid_indices = flux_density_jy >= flux_limit
    ra = ra[valid_indices]
    dec = dec[valid_indices]
    flux_density_jy = flux_density_jy[valid_indices]

    # Assign a default spectral index for GSM2008 sources
    spectral_index = 0.0

    # Create source dictionaries
    sources = [
        {
            "coords": SkyCoord(ra=r * u.deg, dec=d * u.deg, frame="icrs"),
            "flux": f,
            "spectral_index": spectral_index,
        }
        for r, d, f in zip(ra, dec, flux_density_jy)
    ]

    print(
        f"GSM2008 sources loaded with {len(sources)} sources after applying flux limit of {flux_limit} Jy."
    )
    
    # Measure the memory size of the GSM2008 sources
    gsm_size = asizeof.asizeof(sources)
    print(
        f"GSM2008 sources loaded with {len(sources)} sources after applying flux limit of {flux_limit} Jy."
        f" Total size: {gsm_size / (1024**2):.2f} MB"
    )
    return sources




def get_sources(use_gleam=False, use_gsm=False, flux_limit=1.0, frequency=76e6, nside=32):
    """
    Returns either test sources, GLEAM catalog sources, or GSM2008 sources based on the arguments.

    Parameters:
    - use_gleam (bool): If True, loads the GLEAM catalog.
    - use_gsm (bool): If True, loads GSM2008 sources.
    - flux_limit (float): Flux limit for filtering GLEAM sources.
    - frequency (float): Frequency for GSM sources in Hz.
    - nside (int): Healpix resolution parameter for GSM sources.

    Returns:
    list: List of sources with "coords", "flux", and "spectral_index".
    """
    if use_gleam:
        print("Loading GLEAM catalog...")
        return load_gleam_catalog(flux_limit)
    elif use_gsm:
        print("Loading GSM2008 sources...")
        return load_gsm2008_sources(frequency=frequency, nside=nside, flux_limit=flux_limit)
    else:
        print("Using test sources...")
        return test_sources

