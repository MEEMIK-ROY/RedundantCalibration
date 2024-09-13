# src/source.py
from astropy.coordinates import SkyCoord
import astropy.units as u

sources = [
    {'coords': SkyCoord(ra=0*u.deg, dec=-30.72152777777791*u.deg), 'flux': 2},
    {'coords': SkyCoord(ra=120*u.deg, dec=-30.72152777777791*u.deg), 'flux': 4},
    {'coords': SkyCoord(ra=240*u.deg, dec=-30.72152777777791*u.deg), 'flux': 6}
]

sp_index = 0  # Spectral index
