# src/main.py

"""
Entry point for the RRI visibility simulator. Handles antenna generation, baseline calculations,
sky model integration, visibility computations, and visualization.
"""

# Standard library imports
import os
import argparse
import time
import json
import numpy as np
from astropy.constants import c
from astropy.time import TimeDelta
import astropy.units as u
import matplotlib.pyplot as plt
from bokeh.io import show
import h5py

# Custom module imports
from antenna import generate_antennas, read_antenna_positions
from baseline import generate_baselines
from source import get_sources
from observation import get_location_and_time
from visibility import (
    calculate_visibility_original,
    calculate_visibility_optimized,
    calculate_polarized_visibility,
    calculate_modulus_phase,
)
from skymodel import collect_sky_model_for_time
from plot import plot_visibility, plot_heatmaps, plot_modulus_vs_frequency


def main():
    """
    Main function for the simulator. Parses arguments, initializes parameters,
    computes visibilities, and generates visualizations.
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Radio Astronomy Visibility Simulation"
    )
    parser.add_argument(
        "--num_antennas", type=int, default=3, help="Number of antennas"
    )
    parser.add_argument(
        "--spacing", type=float, default=14.0, help="Spacing between antennas in meters"
    )
    parser.add_argument(
        "--lat", type=float, help="Latitude of the observation location in degrees"
    )
    parser.add_argument(
        "--lon", type=float, help="Longitude of the observation location in degrees"
    )
    parser.add_argument(
        "--height", type=float, help="Height of the observation location in meters"
    )
    parser.add_argument(
        "--antenna_positions_file",
        type=str,
        help="File containing antenna positions. Each line should contain x y z coordinates of an antenna, separated by spaces.",
    )
    parser.add_argument(
        "--use_non_optimized",
        action="store_true",
        help="Use the original non-optimized visibility calculation.",
    )
    parser.add_argument(
        "--use_gleam",
        action="store_true",
        help="Use the GLEAM catalog as the source model.",
    )
    parser.add_argument(
        "--flux_limit",
        type=float,
        default=1.0,
        help="Flux limit for the GLEAM sources (Jy).",
    )
    parser.add_argument(
        "--gleam_flux_limit",
        type=float,
        default=1.0,
        help="Flux limit for GLEAM sources.",
    )
    parser.add_argument(
        "--plotting",
        type=str,
        choices=["matplotlib", "bokeh"],
        default="bokeh",
        help='Choose the plotting library: "matplotlib" or "bokeh".',
    )
    parser.add_argument(
        "--use_polarization",
        action="store_true",
        help="Include polarization in the visibility calculation.",
    )
    parser.add_argument(
        "--plot_skymodel_every_hour",
        action="store_true",
        help="Plot the sky model at each hour during the simulation.",
    )
    parser.add_argument(
        "--skymodel_frequency",
        type=float,
        default=150.0,
        help="Frequency in MHz for the sky model.",
    )
    parser.add_argument(
        "--save_skymodel_plots",
        action="store_true",
        help="Save the sky model plots to files.",
    )
    parser.add_argument(
        "--skymodel_plot_dir",
        type=str,
        default="skymodel_plots",
        help="Directory to save sky model plots.",
    )
    parser.add_argument(
        "--fov_radius_deg",
        type=float,
        default=5.0,
        help="Field of view radius in degrees.",
    )
    parser.add_argument(
        "--theta_HPBW",
        type=float,
        default=10.9,
        help="Half Power Beam Width (HPBW) of the antenna in degrees.",
    )
    parser.add_argument(
        "--time_interval",
        type=float,
        default=1.0,
        help="Time interval between calculations.",
    )
    parser.add_argument(
        "--time_interval_unit",
        type=str,
        choices=["hours", "seconds"],
        default="hours",
        help='Unit for time interval: "hours" or "seconds".',
    )
    parser.add_argument(
        "--total_duration",
        type=float,
        default=24.0,
        help="Total duration of the calculation.",
    )
    parser.add_argument(
        "--total_duration_unit",
        type=str,
        choices=["days", "hours", "seconds"],
        default="hours",
        help='Unit for total duration: "days", "hours", or "seconds".',
    )
    parser.add_argument(
        "--output_file",
        type=str,
        default="visibility_data.h5",
        help="Output HDF5 file to store visibility data.",
    )

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
        antennas = generate_antennas(
            num_antennas=args.num_antennas, spacing=args.spacing
        )

    # Generate baselines between antennas
    baselines = generate_baselines(antennas)

    # Get observation location and start time
    location, obstime_start = get_location_and_time(args.lat, args.lon, args.height)

    # Lists to store sky model plotting functions and figures
    sky_model_plots = []
    sky_model_figures = []

    # Create directory for sky model plots if saving is enabled
    if args.save_skymodel_plots:
        os.makedirs(args.skymodel_plot_dir, exist_ok=True)

    # Load sources (either test sources or GLEAM catalog)
    sources = get_sources(use_gleam=args.use_gleam, flux_limit=args.gleam_flux_limit)

    # Frequency and wavelength setup
    freqs = np.linspace(1e8, 2e8, 100)  # Frequencies from 100 MHz to 200 MHz
    wavelengths = c / freqs

    # Convert theta_HPBW from degrees to radians
    theta_HPBW = np.deg2rad(args.theta_HPBW)

    # Convert time_interval to seconds
    if args.time_interval_unit == "hours":
        time_interval_seconds = args.time_interval * 3600
    elif args.time_interval_unit == "seconds":
        time_interval_seconds = args.time_interval
    else:
        raise ValueError("Invalid time_interval_unit.")

    # Convert total_duration to seconds
    if args.total_duration_unit == "days":
        total_duration_seconds = args.total_duration * 86400
    elif args.total_duration_unit == "hours":
        total_duration_seconds = args.total_duration * 3600
    elif args.total_duration_unit == "seconds":
        total_duration_seconds = args.total_duration
    else:
        raise ValueError("Invalid total_duration_unit.")

    # Time points for simulation
    time_points = np.arange(0, total_duration_seconds, time_interval_seconds)

    # Convert time points to MJD format
    mjd_time_points = (obstime_start + TimeDelta(time_points, format="sec")).mjd

    # Initialize dictionaries to store modulus and phase over time for each baseline
    moduli_over_time = {key: [] for key in baselines.keys()}
    phases_over_time = {key: [] for key in baselines.keys()}

    start_time = time.time()

    # # Filename to save visibility results
    # visibility_filename = "visibility_results.json"

    # # Clear or create the visibility file
    # with open(visibility_filename, "w") as f:
    #     f.write("")  # Empty the file

    # Choose the visibility calculation method based on arguments
    if args.use_polarization:
        print("Using polarized visibility calculation.")
        calculate_visibility = calculate_polarized_visibility
    else:
        if args.use_non_optimized:
            print("Using original non-optimized visibility calculation.")
            calculate_visibility = calculate_visibility_original
        else:
            print("Using optimized non-polarized visibility calculation.")
            calculate_visibility = calculate_visibility_optimized

    # Initialize dictionaries to store complex visibility over time for each baseline
    visibilities = {key: [] for key in baselines.keys()}

    
    # Loop over time points to calculate visibility
    for idx, current_time in enumerate(time_points):
        # Update observation time
        obstime = obstime_start + TimeDelta(current_time, format="sec")

        # Calculate visibility for current time
        visibility_dict = calculate_visibility(
            antennas,
            baselines,
            sources,
            location,
            obstime,
            wavelengths,
            freqs,
            theta_HPBW,
        )

        # Append visibility data for each baseline
        for key in baselines.keys():
            visibilities[key].append(visibility_dict[key])

        # # Collect sky model plotting functions if required
        if args.plot_skymodel_every_hour:
            plot_function = collect_sky_model_for_time(
                location=location,
                obstime=obstime,
                frequency=args.skymodel_frequency,
                fov_radius_deg=args.fov_radius_deg,
                gleam_sources=sources if args.use_gleam else None,
            )
            fig = plot_function()
            plt.close(fig)
            sky_model_figures.append(fig)

        # # Calculate modulus and phase of visibility
        moduli, phases = calculate_modulus_phase(visibility_dict)

        # # Append modulus and phase to the time series data
        for key in baselines.keys():
            moduli_over_time[key].append(moduli[key])
            phases_over_time[key].append(phases[key])

        # Save visibility data to file
        # visibility_dict_to_save = {
        #     str(key): {"modulus": moduli[key].tolist(), "phase": phases[key].tolist()}
        #     for key in visibility_dict
        # }

        # with open(visibility_filename, "a") as f:
        #     json.dump({"time": current_time, "visibility": visibility_dict_to_save}, f)
        #     f.write("\n")  # Newline for each time point's data

        # Progress reporting every 50 steps
        if idx % 200 == 0:
            elapsed_time = time.time() - start_time
            if idx > 0:
                time_per_step = elapsed_time / idx
                remaining_steps = len(time_points) - idx
                estimated_remaining_time = time_per_step * remaining_steps
                print(
                    f"Time step {idx+1}/{len(time_points)}: "
                    f"Estimated remaining time: {estimated_remaining_time:.2f} seconds"
                )
            else:
                print(
                    f"Time step {idx+1}/{len(time_points)}: Time elapsed: {elapsed_time:.2f} seconds"
                )

    #    # **Modified code to save data to HDF5**
    #     with h5py.File(args.output_file, "w") as h5file:
    #         # Save complex visibility as primary dataset
    #         for key, vis in visibilities.items():
    #             # Stack the list of arrays into a 2D NumPy array
    #             vis_array = np.stack(vis)  # Shape: (number_of_time_points, number_of_frequencies)
    #             # Ensure the data type is complex128
    #             vis_array = vis_array.astype(np.complex128)
    #             print(f"Key: {key}, vis_array.shape: {vis_array.shape}, dtype: {vis_array.dtype}")
    #             dset_name = f"visibility_{key[0]}_{key[1]}"
    #             h5file.create_dataset(dset_name, data=vis_array, dtype="complex128")

    #         # Save metadata as datasets or attributes
    #         h5file.create_dataset("frequencies", data=freqs)
    #         h5file.create_dataset("time_points", data=time_points)
    #         h5file.create_dataset("wavelengths", data=wavelengths.value)
    #         # Compute LST values for each time point
    #         lst_array = (obstime_start + TimeDelta(time_points, format="sec")).sidereal_time("apparent", location.lon).hour
    #         h5file.create_dataset("lst", data=lst_array)
    #         h5file.attrs["flux_limit"] = args.flux_limit
    #         h5file.attrs["gleam_flux_limit"] = args.gleam_flux_limit
    #         h5file.attrs["num_antennas"] = args.num_antennas
    #         h5file.attrs["spacing"] = args.spacing
    #         h5file.attrs["location"] = str(location)
    #         h5file.attrs["theta_HPBW"] = args.theta_HPBW
    # ## also save the number of point sources
    # print(f"Simulation data saved to {args.output_file}")

    # Save computed data to an HDF5 file if output_file argument is provided
    if args.output_file:
        with h5py.File(args.output_file, "w") as h5file:
            # Create a group for each baseline
            for key, vis in visibilities.items():
                baseline_group = h5file.create_group(
                    f"baseline_{key}"
                )  # Create a group for the baseline

                # Save complex visibility
                vis_array = np.stack(vis)  # Convert to 2D NumPy array
                baseline_group.create_dataset(
                    "complex_visibility",
                    data=vis_array.astype(np.complex128),
                    dtype="complex128",
                )

            # Save frequencies and time points
            h5file.create_dataset("frequencies", data=freqs)
            h5file.create_dataset("time_points_mjd", data=mjd_time_points)

            # Save metadata as attributes
            h5file.attrs["gleam_flux_limit"] = args.gleam_flux_limit
            h5file.attrs["num_antennas"] = len(
                antennas
            )  # Calculated from the antennas list
            h5file.attrs["num_baselines"] = len(
                baselines
            )  # Calculated from the baselines dictionary
            h5file.attrs["location"] = str(location)
            h5file.attrs["theta_HPBW"] = args.theta_HPBW

        print(f"Simulation data saved to {args.output_file}")

    # Calculate total memory usage for visibilities in MB
    total_memory_mb = 0
    for key in visibilities.keys():
        # Convert the list of arrays for each baseline into a stacked 2D NumPy array
        vis_array = np.stack(visibilities[key])  # Shape: (time_steps, frequencies)
        total_memory_mb += vis_array.nbytes / (1024**2)  # Convert bytes to MB

    print(f"Total memory used by visibility data: {total_memory_mb:.2f} MB")

    # After the loop, execute all the sky model plot functions
    # if args.plot_skymodel_every_hour:
    #     for i, plot_func in enumerate(sky_model_plots):
    #         plt.figure()  # Create a new figure for each plot
    #         plot_func()  # Generate each sky model plot

    # Convert lists to numpy arrays for plotting
    for key in baselines.keys():
        moduli_over_time[key] = np.array(moduli_over_time[key])
        phases_over_time[key] = np.array(phases_over_time[key])

    # Convert time_points to desired units for plotting, e.g., hours
    # time_points_hours = time_points / 3600  # Convert seconds to hours for plotting

    # Generate plots based on the selected plotting library
    fig1 = plot_visibility(
        moduli_over_time,
        phases_over_time,
        baselines,
        mjd_time_points,  # Use time in hours for plotting
        freqs,
        total_duration_seconds,
        plotting=args.plotting,
    )
    fig2 = plot_heatmaps(
        moduli_over_time,
        phases_over_time,
        baselines,
        freqs,
        total_duration_seconds,
        mjd_time_points,
        plotting=args.plotting,
    )
    fig3 = plot_modulus_vs_frequency(
        moduli_over_time, phases_over_time, baselines, freqs, 0, plotting=args.plotting
    )

    # Display plots
    if args.plotting == "bokeh":
        show(fig1)
        show(fig3)
        show(fig2)
    else:
        plt.show()


if __name__ == "__main__":
    main()
