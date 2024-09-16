# src/main.py

import numpy as np
from astropy.constants import c
from antenna import generate_antennas, read_antenna_positions
from baseline import generate_baselines
from source import get_sources
from observation import get_location_and_time
from visibility import (
    calculate_visibility_original,
    calculate_visibility_optimized,
    calculate_polarized_visibility,
    calculate_modulus_phase
)
from plot import plot_visibility, plot_heatmaps, plot_modulus_vs_frequency
from astropy.time import TimeDelta
import argparse
import astropy.units as u
import matplotlib.pyplot as plt
import time
import json
from bokeh.io import show


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Radio Astronomy Visibility Simulation')
    parser.add_argument('--num_antennas', type=int, default=3, help='Number of antennas')
    parser.add_argument('--spacing', type=float, default=14.0, help='Spacing between antennas in meters')
    parser.add_argument('--lat', type=float, help='Latitude of the observation location in degrees')
    parser.add_argument('--lon', type=float, help='Longitude of the observation location in degrees')
    parser.add_argument('--height', type=float, help='Height of the observation location in meters')
    parser.add_argument('--sp_index', type=float, default=0.0, help='Spectral index for source flux adjustment.')
    parser.add_argument('--antenna_positions_file', type=str, help='File containing antenna positions. Each line should contain x y z coordinates of an antenna, separated by spaces.')
    parser.add_argument('--use_optimized', action='store_true', help='Use the optimized visibility calculation.')
    parser.add_argument('--use_gleam', action='store_true', help='Use the GLEAM catalog as the source model.')
    parser.add_argument('--flux_limit', type=float, default=1.0, help='Flux limit for the GLEAM sources (Jy).')
    parser.add_argument('--plotting', type=str, choices=['matplotlib', 'bokeh'], default='matplotlib', help='Choose the plotting library: "matplotlib" or "bokeh".')
    parser.add_argument('--use_polarization', action='store_true', help='Include polarization in the visibility calculation.')


    args = parser.parse_args()

    # Input validation
    if args.num_antennas < 2:
        parser.error("Number of antennas must be at least 2.")

    if args.spacing <= 0:
        parser.error("Spacing between antennas must be a positive number.")


    # Generate antennas
    if args.antenna_positions_file:
        # Read antenna positions from file
        antennas = read_antenna_positions(args.antenna_positions_file)
    else:
        # Generate antennas with specified spacing
        antennas = generate_antennas(num_antennas=args.num_antennas, spacing=args.spacing)


    # Generate baselines
    baselines = generate_baselines(antennas)

    # Get location and observation time
    location, obstime_start = get_location_and_time(args.lat, args.lon, args.height)

    # Load sources (either test or GLEAM catalog)
    sources = get_sources(use_gleam=args.use_gleam, flux_limit=args.flux_limit)


    # Frequency and wavelength setup
    freqs = np.linspace(1e8, 2e8, 100)  # Frequencies from 100 MHz to 200 MHz
    wavelengths = c / freqs

    # Time settings
    hours_per_day = 24
    time_interval = 1  # Run every hour
    total_hours = hours_per_day * time_interval

    moduli_over_time = {key: [] for key in baselines.keys()}
    phases_over_time = {key: [] for key in baselines.keys()}

    start_time = time.time()

    # Visibility results filename
    visibility_filename = "visibility_results.json"

    # Clear or create the visibility file
    with open(visibility_filename, 'w') as f:
        f.write("")  # Empty the file

    # Choose the visibility calculation method
    if args.use_polarization:
        print("Using polarized visibility calculation.")
        calculate_visibility = calculate_polarized_visibility
    else:
        if args.use_optimized:
            print("Using optimized non-polarized visibility calculation.")
            calculate_visibility = calculate_visibility_optimized
        else:
            print("Using original non-polarized visibility calculation.")
            calculate_visibility = calculate_visibility_original



    # Calculate visibility over time (every hour)
    for hour in range(0, total_hours, time_interval):
        obstime = obstime_start + TimeDelta(hour * 3600, format='sec')  # Convert hours to seconds
        visibility_dict = calculate_visibility(
            antennas, baselines, sources, location, obstime, wavelengths, freqs, args.sp_index
        )

        # # Calculate and log the elapsed time
        # elapsed_time = time.time() - start_time
        # remaining_hours = total_hours - hour - 1
        # estimated_total_time = elapsed_time / (hour + 1) * total_hours
        # remaining_time = estimated_total_time - elapsed_time
        
        # print(f"Hour {hour+1}/{total_hours}: Estimated remaining time: {remaining_time:.2f} seconds")

        moduli, phases = calculate_modulus_phase(visibility_dict)

        for key in baselines.keys():
            moduli_over_time[key].append(moduli[key])
            phases_over_time[key].append(phases[key])


        # Save visibility continuously
        visibility_dict_to_save = {
            str(key): {
                'modulus': moduli[key].tolist(),
                'phase': phases[key].tolist()
            } for key in visibility_dict
        }

        with open(visibility_filename, 'a') as f:
            json.dump({'hour': hour, 'visibility': visibility_dict_to_save}, f)
            f.write("\n")  # Write newline for each hour's data

        # Timing and estimation
        elapsed_time = time.time() - start_time
        if hour > 0:
            time_per_hour = elapsed_time / hour
            remaining_hours = total_hours - hour
            estimated_remaining_time = time_per_hour * remaining_hours
            print(f"Hour {hour+1}/{total_hours}: Estimated remaining time: {estimated_remaining_time:.2f} seconds")
        else:
            print(f"Hour {hour+1}/{total_hours}: Time elapsed: {elapsed_time:.2f} seconds")

    # Convert lists to numpy arrays
    for key in baselines.keys():
        moduli_over_time[key] = np.array(moduli_over_time[key])
        phases_over_time[key] = np.array(phases_over_time[key])

    # Time points for plotting (in hours)
    time_points = np.arange(0, total_hours, time_interval)

    # Plot based on the selected plotting library
    fig1 = plot_visibility(moduli_over_time, phases_over_time, baselines, time_points, freqs, total_hours, plotting=args.plotting)
    fig2 = plot_heatmaps(moduli_over_time, phases_over_time, baselines, freqs, total_hours, plotting=args.plotting)
    fig3 = plot_modulus_vs_frequency(moduli_over_time, baselines, freqs, 0, plotting=args.plotting)

    if args.plotting == 'bokeh':
        show(fig1)
        show(fig2)
        show(fig3)
    else:
        plt.show()

if __name__ == '__main__':
    main()
