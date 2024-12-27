# src/observation.py

from astropy.coordinates import EarthLocation
from astropy.time import Time
import astropy.units as u

def get_location_and_time(lat=None, lon=None, height=None, starttime=None):
    """
    Returns the observation location and start time.

    Parameters:
    lat (float): Latitude in degrees. Defaults to HERA latitude.
    lon (float): Longitude in degrees. Defaults to HERA longitude.
    height (float): Height in meters. Defaults to HERA height.
    starttime (str): Observation start time in ISO format. Defaults to the current UTC time at HERA.


    Returns:
    tuple: EarthLocation object and Time object for observation start time.
    """
    # Default HERA coordinates if not provided
    default_lat = -30.72152777777791
    default_lon = 21.428305555555557
    default_height = 1073.0

    lat = lat if lat is not None else default_lat
    lon = lon if lon is not None else default_lon
    height = height if height is not None else default_height

    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=height * u.m)
    # Default start time to the current UTC time at HERA's location
    if starttime is None:
        obstime_start = Time.now()
    else:
        obstime_start = Time(starttime, format="isot", scale="utc")

    return location, obstime_start
