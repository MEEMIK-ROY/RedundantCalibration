# src/observation.py

from astropy.coordinates import EarthLocation
from astropy.time import Time

hera_lon = 21.428305555555557
hera_lat = -30.72152777777791
hera_height = 1073.0

location = EarthLocation.from_geodetic(lat=hera_lat, lon=hera_lon, height=hera_height)
obstime_start = Time("2024-06-05T00:00:00", format="isot", scale="utc")
