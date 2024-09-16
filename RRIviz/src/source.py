# src/source.py

from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from pympler import asizeof

# Test sources (used when not loading GLEAM catalog)
test_sources = [
    {'coords': SkyCoord(ra=0*u.deg, dec=-30.72152777777791*u.deg), 'flux': 2},
    {'coords': SkyCoord(ra=120*u.deg, dec=-30.72152777777791*u.deg), 'flux': 4},
    {'coords': SkyCoord(ra=240*u.deg, dec=-30.72152777777791*u.deg), 'flux': 6}
]

sp_index = 0  # Spectral index for both models

# GLEAM catalog loading function
def load_gleam_catalog(flux_limit=1.0):
    """
    Loads the GLEAM Extragalactic catalog from VizieR using astroquery.
    Applies a flux limit to reduce the number of sources.

    Parameters:
    flux_limit (float): Minimum flux (in Jy) to include a source.

    Returns:
    list: List of sources with their coordinates and fluxes.
    """
    # Set up Vizier query
    Vizier.ROW_LIMIT = -1  # No row limit
    # Retrieve the GLEAM Extragalactic Catalog
    catalog_list = Vizier.find_catalogs('VIII/100/gleamegc')
    if not catalog_list:
        raise Exception("GLEAM catalog not found in VizieR.")
    catalog = Vizier.get_catalogs('VIII/100/gleamegc')[0]

    # Apply flux limit and extract positions and fluxes
    sources = []
    for row in catalog:
        flux = row['Fpwide']  # Wide-band peak flux density in Jy
        if flux >= flux_limit:
            ra = row['RAJ2000'] * u.deg
            dec = row['DEJ2000'] * u.deg
            coords = SkyCoord(ra=ra, dec=dec)
            sources.append({'coords': coords, 'flux': flux})



    # Measure the memory size of the loaded GLEAM catalog
    gleam_size = asizeof.asizeof(sources)
    print(f"GLEAM catalog loaded with {len(sources)} sources. Total size: {gleam_size / (1024**2):.2f} MB")

    return sources

def get_sources(use_gleam=False, flux_limit=1.0):
    """
    Returns either the test sources or GLEAM catalog sources based on the argument.

    Parameters:
    use_gleam (bool): If True, loads the GLEAM catalog. Otherwise, uses the test sources.
    flux_limit (float): Flux limit for filtering GLEAM sources.

    Returns:
    list: List of sources with their coordinates and fluxes.
    """
    if use_gleam:
        print("Loading GLEAM catalog...")
        return load_gleam_catalog(flux_limit)
    else:
        print("Using test sources...")
        return test_sources
