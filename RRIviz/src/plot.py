# src/plot.py

from bokeh.plotting import figure, show
from bokeh.palettes import Viridis256, Turbo256
from bokeh.layouts import column
from bokeh.models import DatetimeTicker, HoverTool, Legend
from bokeh.io import output_file, save, reset_output
import matplotlib.pyplot as plt
import numpy as np
import os
from bokeh.resources import CDN
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta


def plot_visibility(
    moduli_over_time,
    phases_over_time,
    baselines,
    mjd_time_points,
    freqs,
    total_seconds,
    plotting="bokeh",
    save_simulation_data=False,
    folder_path=None,
):
    """
    Create plots of the modulus and phase of visibility versus time for each baseline.

    Parameters:
    moduli_over_time (dict): Dictionary of visibility moduli over time for each baseline.
    phases_over_time (dict): Dictionary of visibility phases over time for each baseline.
    baselines (dict): Dictionary of baselines between antennas.
    time_points (ndarray): Array of time points corresponding to the observations.
    freqs (ndarray): Array of frequencies used in the observations.
    total_seconds (float): Total duration of the observation in seconds.
    plotting (str): Plotting library to use ('matplotlib' or 'bokeh').

    Returns:
    Figure object: The generated plot figure(s) based on the specified plotting library.
    """
    
    # Convert MJD to human-readable datetime
    time_points_datetime = Time(mjd_time_points, format="mjd").to_datetime()

    baseline_keys = list(baselines.keys())
    colors = Turbo256

    if plotting == "bokeh":
        # Bokeh plotting
        plots = []
        for key in baseline_keys:
            # # Check if polarization is included by examining the shape of the data
            # if moduli_over_time[key].shape[-1] == 2:
            #     # Ex component
            #     p_mod_ex = figure(
            #         width=800,
            #         height=300,
            #         title=f"Modulus of Visibility (Ex) vs Time for Baseline {key}",
            #     )
            #     p_mod_ex.line(
            #         time_points,
            #         moduli_over_time[key][:, 0, 0],
            #         line_width=2,
            #         legend_label=f"Ex Baseline {key}",
            #     )
            #     p_mod_ex.xaxis.axis_label = "Time (seconds)"
            #     p_mod_ex.yaxis.axis_label = "Modulus of Visibility"
            #     p_mod_ex.legend.location = "top_left"

            #     # Ey component
            #     p_mod_ey = figure(
            #         width=800,
            #         height=300,
            #         title=f"Modulus of Visibility (Ey) vs Time for Baseline {key}",
            #     )
            #     p_mod_ey.line(
            #         time_points,
            #         moduli_over_time[key][:, 0, 1],
            #         line_width=2,
            #         legend_label=f"Ey Baseline {key}",
            #     )
            #     p_mod_ey.xaxis.axis_label = "Time (seconds)"
            #     p_mod_ey.yaxis.axis_label = "Modulus of Visibility"
            #     p_mod_ey.legend.location = "top_left"

            #     # Phase components for Ex and Ey
            #     p_phase_ex = figure(
            #         width=800,
            #         height=300,
            #         title=f"Phase of Visibility (Ex) vs Time for Baseline {key}",
            #     )
            #     p_phase_ex.line(
            #         time_points,
            #         phases_over_time[key][:, 0, 0],
            #         line_width=2,
            #         legend_label=f"Ex Baseline {key}",
            #     )
            #     p_phase_ex.xaxis.axis_label = "Time (seconds)"
            #     p_phase_ex.yaxis.axis_label = "Phase of Visibility (radians)"
            #     p_phase_ex.legend.location = "top_left"

            #     p_phase_ey = figure(
            #         width=800,
            #         height=300,
            #         title=f"Phase of Visibility (Ey) vs Time for Baseline {key}",
            #     )
            #     p_phase_ey.line(
            #         time_points,
            #         phases_over_time[key][:, 0, 1],
            #         line_width=2,
            #         legend_label=f"Ey Baseline {key}",
            #     )
            #     p_phase_ey.xaxis.axis_label = "Time (seconds)"
            #     p_phase_ey.yaxis.axis_label = "Phase of Visibility (radians)"
            #     p_phase_ey.legend.location = "top_left"

            #     plots.extend([p_mod_ex, p_mod_ey, p_phase_ex, p_phase_ey])
            # else:
            
            
            # Define a ticker with finer granularity (e.g., hourly ticks)
            datetime_ticker = DatetimeTicker(desired_num_ticks=12)  # Adjust `desired_num_ticks` as needed

            # No polarization, single component
            p_mod = figure(
                width=800,
                height=300,
                title=f"Modulus of Visibility vs Time for Baseline {key}",
                x_axis_type="datetime",
            )
            p_mod.line(
                time_points_datetime,
                moduli_over_time[key][:, 0],
                line_width=2,
                legend_label=f"Baseline {key}",
            )
            p_mod.xaxis.axis_label = "Time"
            p_mod.yaxis.axis_label = "Modulus of Visibility"
            p_mod.legend.location = "top_left"
            p_mod.xaxis.ticker = datetime_ticker


            p_phase = figure(
                width=800,
                height=300,
                title=f"Phase of Visibility vs Time for Baseline {key}",
                x_axis_type="datetime",
            )
            p_phase.line(
                time_points_datetime,
                np.unwrap(phases_over_time[key][:, 0]),
                line_width=2,
                legend_label=f"Baseline {key}",
            )
            p_phase.xaxis.axis_label = "Time"
            p_phase.yaxis.axis_label = "Phase of Visibility (radians)"
            p_phase.legend.location = "top_left"
            p_phase.xaxis.ticker = datetime_ticker


            plots.append(p_mod)
            plots.append(p_phase)
            
        combined_mod = figure(
            width=1400,
            height=1400,
            title="Modulus of Visibility vs Time for All Baselines",
            x_axis_type="datetime",
        )

        lines = []
        for idx, key in enumerate(baseline_keys):
            color = colors[int((idx / len(baseline_keys)) * 255)]
            line = combined_mod.line(
                time_points_datetime,
                moduli_over_time[key][:, 0],
                line_width=2,
                color=color,
                name=str(key),
            )
            lines.append((f"Baseline {key}", [line]))

            hover = HoverTool(
                renderers=[line],
                tooltips=[
                    ("Time", "@x{%F %T}"),
                    ("Value", "@y"),
                    ("Baseline", str(key)),
                ],
                formatters={"@x": "datetime"},
                mode="mouse",
            )
            combined_mod.add_tools(hover)

        combined_mod.xaxis.axis_label = "Time"
        combined_mod.yaxis.axis_label = "Modulus of Visibility"
        combined_mod.xaxis.ticker = DatetimeTicker(desired_num_ticks=12)

        legend = Legend(
            items=lines,
            location="center",
            click_policy="hide",
            title="Baselines",
        )
        legend.ncols = 10
        combined_mod.add_layout(legend, "below")

        # combined_mod.xaxis.axis_label = "Time"
        # combined_mod.yaxis.axis_label = "Modulus of Visibility"
        # combined_mod.legend.label_text_font_size = "8pt"
        # combined_mod.legend.spacing = 1
        # combined_mod.legend.location = "top_left"
        # combined_mod.legend.click_policy = "hide"  # Allow toggling baselines on/off
        # combined_mod.xaxis.ticker = DatetimeTicker(desired_num_ticks=12)


        # Combined Phase vs Time
        combined_phase = figure(
            width=1400,
            height=1400,
            title="Combined Phase of Visibility vs Time for All Baselines",
            x_axis_type="datetime",
        )
        lines = []
        for idx, key in enumerate(baseline_keys):
            color = colors[int((idx / len(baseline_keys)) * 255)]
            line = combined_phase.line(
                time_points_datetime,
                np.unwrap(phases_over_time[key][:, 0]),
                line_width=2,
                color=color,
                # legend_label=f"Baseline {key}",
                name=str(key),  # Convert key to string
            )
            lines.append((f"Baseline {key}", [line]))

            # Add hover tool for this line
            hover = HoverTool(
                renderers=[line],
                tooltips=[
                    ("Time", "@x{%F %T}"),  # Time in human-readable format
                    ("Value", "@y"),        # Value (phase)
                    ("Baseline", str(key)), # Convert key to string for tooltip
                ],
                formatters={"@x": "datetime"},  # Formatter for time
                mode="mouse",  
            )
            combined_phase.add_tools(hover)

        combined_phase.xaxis.axis_label = "Time"
        combined_phase.yaxis.axis_label = "Phase of Visibility (radians)"
        # combined_phase.legend.location = "top_left"
        # combined_phase.legend.label_text_font_size = "8pt"
        # combined_phase.legend.spacing = 1
        # combined_phase.legend.click_policy = "hide"  # Allow toggling baselines on/off
        combined_phase.xaxis.ticker = DatetimeTicker(desired_num_ticks=12)

        # Create a scrollable legend
        legend = Legend(
            items=lines,  # Use the collected line renderers
            location="center",
            click_policy="hide",  # Allow clicking to hide/show lines
            title="Baselines",
        )
        
        legend.ncols = 10
        
        # Restrict legend height to make it scrollable
        combined_phase.add_layout(legend, "below")


        plots.append(combined_mod)
        plots.append(combined_phase)
        

        # Combine plots into a column
        plot_column = column(*plots)

        # Save the entire column if required
        if save_simulation_data and folder_path:
            file_path = os.path.join(folder_path, "visibility-phase-lsts.html")
            save(plot_column, filename=file_path, resources=CDN, title="Visibility/Phase Plots")
            print(f"Saved visibility plots column to {file_path}")

        return plot_column

    else:
        return None
    #     # Matplotlib plotting
    #     num_baselines = len(baselines)
    #     fig1, ax1 = plt.subplots(num_baselines, 2, figsize=(15, 5 * num_baselines))

    #     for i, key in enumerate(baseline_keys):
    #         # # If polarization is included, plot both Ex and Ey components
    #         # if moduli_over_time[key].shape[-1] == 2:
    #         #     ax1[i, 0].plot(
    #         #         time_points,
    #         #         moduli_over_time[key][:, 0, 0],
    #         #         marker="o",
    #         #         label=f"Ex Baseline {key}",
    #         #     )
    #         #     ax1[i, 0].plot(
    #         #         time_points,
    #         #         moduli_over_time[key][:, 0, 1],
    #         #         marker="o",
    #         #         label=f"Ey Baseline {key}",
    #         #     )
    #         #     ax1[i, 0].set_xlabel("Time (seconds)")
    #         #     ax1[i, 0].set_ylabel("Modulus of Visibility")
    #         #     ax1[i, 0].set_title(
    #         #         f"Modulus of Visibility (Ex, Ey) vs Time for Baseline {key}"
    #         #     )
    #         #     ax1[i, 0].legend()
    #         #     ax1[i, 0].grid(True)

    #         #     ax1[i, 1].plot(
    #         #         time_points,
    #         #         phases_over_time[key][:, 0, 0],
    #         #         marker="o",
    #         #         label=f"Ex Baseline {key}",
    #         #     )
    #         #     ax1[i, 1].plot(
    #         #         time_points,
    #         #         phases_over_time[key][:, 0, 1],
    #         #         marker="o",
    #         #         label=f"Ey Baseline {key}",
    #         #     )
    #         #     ax1[i, 1].set_xlabel("Time (seconds)")
    #         #     ax1[i, 1].set_ylabel("Phase of Visibility (radians)")
    #         #     ax1[i, 1].set_title(
    #         #         f"Phase of Visibility (Ex, Ey) vs Time for Baseline {key}"
    #         #     )
    #         #     ax1[i, 1].legend()
    #         #     ax1[i, 1].grid(True)
    #         # else:
    #         # No polarization, single component
    #         # ax1[i, 0].plot(
    #         #     time_points,
    #         #     moduli_over_time[key][:, 0],
    #         #     marker="o",
    #         #     label=f"Baseline {key}",
    #         # )
    #         # ax1[i, 0].set_xlabel("Time")
    #         # ax1[i, 0].set_ylabel("Modulus of Visibility")
    #         # ax1[i, 0].set_title(f"Modulus of Visibility vs Time for Baseline {key}")
    #         # ax1[i, 0].legend()
    #         # ax1[i, 0].grid(True)

    #         # ax1[i, 1].plot(
    #         #     time_points,
    #         #     np.unwrap(phases_over_time[key][:, 0]),
    #         #     marker="o",
    #         #     label=f"Baseline {key}",
    #         # )
    #         # ax1[i, 1].set_xlabel("Time")
    #         # ax1[i, 1].set_ylabel("Phase of Visibility (radians)")
    #         # ax1[i, 1].set_title(f"Phase of Visibility vs Time for Baseline {key}")
    #         # ax1[i, 1].legend()
    #         # ax1[i, 1].grid(True)

    #     plt.tight_layout()
    #     return fig1


def plot_heatmaps(
    moduli_over_time,
    phases_over_time,
    baselines,
    freqs,
    total_seconds,
    mjd_time_points, 
    plotting="bokeh",
    save_simulation_data=False,
    folder_path=None,
):
    """
    Create heatmaps for visibility modulus and phase over time and frequency.

    Parameters:
    moduli_over_time (dict): Dictionary of visibility moduli over time for each baseline.
    phases_over_time (dict): Dictionary of visibility phases over time for each baseline.
    baselines (dict): Dictionary of baselines between antennas.
    freqs (ndarray): Array of frequencies used in the observations.
    total_seconds (float): Total duration of the observation in seconds.
    plotting (str): Plotting library to use ('matplotlib' or 'bokeh').

    Returns:
    Figure object: The generated heatmap figure(s) based on the specified plotting library.
    """
    
    # Convert MJD to human-readable datetime
    time_points_datetime = Time(mjd_time_points, format="mjd").to_datetime()
    baseline_keys = list(baselines.keys())

    if plotting == "bokeh":
        plots = []
        for key in baseline_keys:
            # if moduli_over_time[key].shape[-1] == 2:
            #     # Combine Ex and Ey components
            #     moduli_total = np.sqrt(
            #         moduli_over_time[key][:, :, 0] ** 2
            #         + moduli_over_time[key][:, :, 1] ** 2
            #     )
            #     phases_total = np.sqrt(
            #         phases_over_time[key][:, :, 0] ** 2
            #         + phases_over_time[key][:, :, 1] ** 2
            #     )
            # else:
            moduli_total = moduli_over_time[key]
            phases_total = np.unwrap(phases_over_time[key], axis=0)
            
            # Ensure data is in the correct format (list of 2D arrays)
            moduli_image = [moduli_total.T]  # Transpose to match Bokeh's image orientation
            phases_image = [phases_total.T]

            # Modulus heatmap
            p_mod = figure(
                width=800,
                height=300,
                title=f"Modulus of Visibility Heatmap for Baseline {key}",
                x_axis_type="datetime",

            )
            p_mod.image(
                image=moduli_image,
                x=time_points_datetime[0],  # Start of the datetime range
                y=freqs[0] / 1e6,  # Start of frequency in MHz
                dw=(time_points_datetime[-1] - time_points_datetime[0]).total_seconds() * 1e3,  # Time duration in ms
                dh=(freqs[-1] - freqs[0]) / 1e6,  # Frequency range in MHz
                palette=Viridis256,
            )
            p_mod.xaxis.axis_label = "Time"
            p_mod.yaxis.axis_label = "Frequency (MHz)"
            p_mod.xaxis.ticker = DatetimeTicker(desired_num_ticks=12)  # Add finer ticks


            # Phase heatmap
            p_phase = figure(
                width=800,
                height=300,
                title=f"Phase of Visibility Heatmap for Baseline {key}",
                x_axis_type="datetime",

            )
            p_phase.image(
                image=phases_image,
                x=time_points_datetime[0],  # Start of the datetime range
                y=freqs[0] / 1e6,  # Start of frequency in MHz
                dw=(time_points_datetime[-1] - time_points_datetime[0]).total_seconds() * 1e3,  # Time duration in ms
                dh=(freqs[-1] - freqs[0]) / 1e6,  # Frequency range in MHz
                palette=Viridis256,
            )
            p_phase.xaxis.axis_label = "Time"
            p_phase.yaxis.axis_label = "Frequency (MHz)"
            p_phase.xaxis.ticker = DatetimeTicker(desired_num_ticks=12)  # Add finer ticks


            plots.append(p_mod)
            plots.append(p_phase)
            
        # Combine plots into a column
        plot_column = column(*plots, sizing_mode="stretch_both")

        # Save the entire column if required
        if save_simulation_data and folder_path:
            file_path = os.path.join(folder_path, "heatmaps-freq-time.html")
            save(plot_column, filename=file_path, resources=CDN, title="Visibility Heatmaps")
            print(f"Saved heatmaps column to {file_path}")

        return plot_column

    else:
        return None
        # # Matplotlib plotting
        # num_baselines = len(baselines)
        # fig2, ax2 = plt.subplots(num_baselines, 2, figsize=(15, 5 * num_baselines))

        # for i, key in enumerate(baseline_keys):
        #     # if moduli_over_time[key].shape[-1] == 2:
        #     #     # Combine Ex and Ey components
        #     #     moduli_total = np.sqrt(
        #     #         moduli_over_time[key][:, :, 0] ** 2
        #     #         + moduli_over_time[key][:, :, 1] ** 2
        #     #     )
        #     #     phases_total = np.sqrt(
        #     #         phases_over_time[key][:, :, 0] ** 2
        #     #         + phases_over_time[key][:, :, 1] ** 2
        #     #     )
        #     # else:
        #     moduli_total = moduli_over_time[key]
        #     phases_total = phases_over_time[key]

        #     # Modulus heatmap
        #     im0 = ax2[i, 0].imshow(
        #         moduli_total,
        #         aspect="auto",
        #         origin="lower",
        #         extent=[freqs[0] / 1e6, freqs[-1] / 1e6, 0, total_seconds],
        #         cmap="twilight",  # Updated color scheme to 'twilight'
        #     )
        #     ax2[i, 0].set_xlabel("Frequency (MHz)")
        #     ax2[i, 0].set_ylabel("Time (seconds)")
        #     ax2[i, 0].set_title(f"Modulus of Visibility for Baseline {key}")
        #     fig2.colorbar(im0, ax=ax2[i, 0])

        #     # Phase heatmap
        #     im1 = ax2[i, 1].imshow(
        #         np.unwrap(phases_total),
        #         aspect="auto",
        #         origin="lower",
        #         extent=[freqs[0] / 1e6, freqs[-1] / 1e6, 0, total_seconds],
        #         cmap="twilight",  # Updated color scheme to 'twilight'
        #     )
        #     ax2[i, 1].set_xlabel("Frequency (MHz)")
        #     ax2[i, 1].set_ylabel("Time (seconds)")
        #     ax2[i, 1].set_title(f"Phase of Visibility for Baseline {key}")
        #     fig2.colorbar(im1, ax=ax2[i, 1])

        # plt.tight_layout()
        # return fig2



def plot_modulus_vs_frequency(
    moduli_over_time, phases_over_time, baselines, freqs, time_index, plotting="bokeh",    save_simulation_data=False,
    folder_path=None,
):
    """
    Create plots of the modulus of visibility versus frequency at a specific time point.

    Parameters:
    moduli_over_time (dict): Dictionary of visibility moduli over time for each baseline.
    phases_over_time (dict): Dictionary of visibility phases over time for each baseline.
    baselines (dict): Dictionary of baselines between antennas.
    freqs (ndarray): Array of frequencies used in the observations.
    time_index (int): Index of the time point at which to plot the modulus.
    plotting (str): Plotting library to use ('matplotlib' or 'bokeh').

    Returns:
    Figure object: The generated plot figure(s) based on the specified plotting library.
    """
    baseline_keys = list(baselines.keys())
    colors = Turbo256

    if plotting == "bokeh":
        plots = []
        for key in baseline_keys:
            
            # Modulus
            # if moduli_over_time[key].shape[-1] == 2:
            #     # Combine Ex and Ey components at the specified time index
            #     moduli_total = np.sqrt(
            #         moduli_over_time[key][time_index, :, 0] ** 2
            #         + moduli_over_time[key][time_index, :, 1] ** 2
            #     )
            # else:
            moduli_total = np.mean(moduli_over_time[key][time_index, :])

            p_mod = figure(
                width=800,
                height=300,
                title=f"Modulus of Visibility vs Frequency for Baseline {key} at Time {time_index} seconds",
            )
            p_mod.line(
                freqs / 1e6, moduli_total, line_width=2, legend_label=f"Baseline {key}"
            )
            p_mod.xaxis.axis_label = "Frequency (MHz)"
            p_mod.yaxis.axis_label = "Modulus of Visibility"
            p_mod.legend.location = "top_left"
            
            # Phase
            phases_total = np.unwrap(phases_over_time[key][time_index, :])
            p_phase = figure(
                width=800,
                height=300,
                title=f"Phase of Visibility vs Frequency for Baseline {key} at Time Index {time_index}",
            )
            p_phase.line(
                freqs / 1e6, phases_total, line_width=2, legend_label=f"Baseline {key}"
            )
            p_phase.xaxis.axis_label = "Frequency (MHz)"
            p_phase.yaxis.axis_label = "Phase of Visibility (radians)"
            p_phase.legend.location = "top_left"

            plots.append(p_mod)
            plots.append(p_phase)
            
        # Combined Modulus vs Frequency 
        combined_mod = figure(
            width=1400,
            height=1400,
            title="Modulus of Visibility vs Frequency for All Baselines",
        )

        lines = []
        for idx, key in enumerate(baseline_keys):
            color = colors[int((idx / len(baseline_keys)) * 255)]
            line = combined_mod.line(
                freqs / 1e6,
                moduli_over_time[key][time_index, :],
                line_width=2,
                color=color,
                name=str(key),
            )
            lines.append((f"Baseline {key}", [line]))

            hover = HoverTool(
                renderers=[line],
                tooltips=[
                    ("Frequency", "@x MHz"),
                    ("Value", "@y"),
                    ("Baseline", str(key)),
                ],
                mode="mouse",
            )
            combined_mod.add_tools(hover)

        combined_mod.xaxis.axis_label = "Frequency (MHz)"
        combined_mod.yaxis.axis_label = "Modulus of Visibility"

        legend = Legend(
            items=lines,
            location="center",
            click_policy="hide",
            title="Baselines",
        )
        legend.ncols = 10
        combined_mod.add_layout(legend, "below")

        # Combined Phase vs Frequency
        combined_phase = figure(
            width=1400,
            height=1400,
            title="Phase of Visibility vs Frequency for All Baselines",
        )

        lines = []
        for idx, key in enumerate(baseline_keys):
            color = colors[int((idx / len(baseline_keys)) * 255)]
            line = combined_phase.line(
                freqs / 1e6,
                np.unwrap(phases_over_time[key][time_index, :]),
                line_width=2,
                color=color,
                name=str(key),
            )
            lines.append((f"Baseline {key}", [line]))

            hover = HoverTool(
                renderers=[line],
                tooltips=[
                    ("Frequency", "@x MHz"),
                    ("Value", "@y radians"),
                    ("Baseline", str(key)),
                ],
                mode="mouse",
            )
            combined_phase.add_tools(hover)

        combined_phase.xaxis.axis_label = "Frequency (MHz)"
        combined_phase.yaxis.axis_label = "Phase of Visibility (radians)"

        legend = Legend(
            items=lines,
            location="center",
            click_policy="hide",
            title="Baselines",
        )
        legend.ncols = 10
        combined_phase.add_layout(legend, "below")

        plots.append(combined_mod)
        plots.append(combined_phase)
        
        plot_column = column(*plots, sizing_mode="stretch_both")

        # Save the entire column if required
        if save_simulation_data and folder_path:
            file_path = os.path.join(folder_path, "modulus-phase-freq.html")
            save(plot_column, filename=file_path, resources=CDN, title="Visibility Modulus/Phase vs Frequency")
            print(f"Saved Modulus vs Frequency column to {file_path}")

        return plot_column

    else:
        return None
        # # Matplotlib plotting
        # num_baselines = len(baselines)
        # fig3, ax3 = plt.subplots(num_baselines, 1, figsize=(15, 5 * num_baselines))

        # for i, key in enumerate(baseline_keys):
        #     # if moduli_over_time[key].shape[-1] == 2:
        #     #     # Combine Ex and Ey components at the specified time index
        #     #     moduli_total = np.sqrt(
        #     #         moduli_over_time[key][time_index, :, 0] ** 2
        #     #         + moduli_over_time[key][time_index, :, 1] ** 2
        #     #     )
        #     # else:
        #     moduli_total = moduli_over_time[key][time_index, :]

        #     ax3[i].plot(freqs / 1e6, moduli_total, marker="o", label=f"Baseline {key}")
        #     ax3[i].set_xlabel("Frequency (MHz)")
        #     ax3[i].set_ylabel("Modulus of Visibility")
        #     ax3[i].set_title(
        #         f"Modulus of Visibility vs Frequency for Baseline {key} at Time {time_index} seconds"
        #     )
        #     ax3[i].legend()
        #     ax3[i].grid(True)

        # plt.tight_layout()
        # return fig3
