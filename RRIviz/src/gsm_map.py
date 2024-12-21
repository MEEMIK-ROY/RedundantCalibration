import pygdsm
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
import astropy.units as u
from datetime import datetime, timedelta

def diffused_sky_model(location, obstime, frequency=76, fov_radius_deg=5, gleam_sources=None):
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
    
    
    return None