# src/main.py

import numpy as np
from astropy.constants import c
from antenna import antennas
from baseline import baselines
from source import sources
from observation import location, obstime_start
from visibility import calculate_visibility, calculate_modulus_phase
from plot import plot_visibility, plot_heatmaps, plot_modulus_vs_frequency
from astropy.time import TimeDelta
import matplotlib.pyplot as plt

# Frequency and wavelength setup
freqs = np.linspace(1e8, 2e8, 100)  # Frequencies from 100 MHz to 200 MHz
wavelengths = c / freqs

# Time settings
hours_per_day = 24
time_interval = 1  # Run every hour
total_hours = hours_per_day * time_interval

moduli_over_time = {key: [] for key in baselines.keys()}
phases_over_time = {key: [] for key in baselines.keys()}

# Calculate visibility over time (every hour)
for hour in range(0, total_hours, time_interval):
    obstime = obstime_start + TimeDelta(hour * 3600, format='sec')  # Convert hours to seconds
    visibility_dict = calculate_visibility(antennas, baselines, sources, location, obstime, wavelengths, freqs, 0)
    moduli, phases = calculate_modulus_phase(visibility_dict)
    
    for key in baselines.keys():
        moduli_over_time[key].append(moduli[key])
        phases_over_time[key].append(phases[key])

# Convert lists to numpy arrays
for key in baselines.keys():
    moduli_over_time[key] = np.array(moduli_over_time[key])
    phases_over_time[key] = np.array(phases_over_time[key])

# Time points for plotting (in hours)
time_points = np.arange(0, total_hours, time_interval)

# Plot all figures at once
fig1 = plot_visibility(moduli_over_time, phases_over_time, baselines, time_points, freqs, total_hours)
fig2 = plot_heatmaps(moduli_over_time, phases_over_time, baselines, freqs, total_hours)
fig3 = plot_modulus_vs_frequency(moduli_over_time, baselines, freqs, 0)

# Display all plots together
plt.show()
