from bokeh.plotting import figure, show
from bokeh.layouts import column
import matplotlib.pyplot as plt
import numpy as np

def plot_visibility(moduli_over_time, phases_over_time, baselines, time_points, freqs, total_seconds, plotting='matplotlib'):
    """
    Create the modulus and phase of visibility vs time using either Matplotlib or Bokeh.
    Returns appropriate plots based on the plotting library specified.
    Handles Ex and Ey polarization components if available.
    """
    baseline_keys = list(baselines.keys())

    if plotting == 'bokeh':
        # Bokeh plotting
        plots = []
        for key in baseline_keys:
            # If polarization is included, plot both Ex and Ey components
            if moduli_over_time[key].shape[-1] == 2:
                # Ex component
                p_mod_ex = figure(width=800, height=300, title=f"Modulus of Visibility (Ex) vs Time for Baseline {key}")
                p_mod_ex.line(time_points, moduli_over_time[key][:, 0, 0], line_width=2, legend_label=f"Ex Baseline {key}")
                p_mod_ex.xaxis.axis_label = 'Time (seconds)'
                p_mod_ex.yaxis.axis_label = 'Modulus of Visibility'
                p_mod_ex.legend.location = "top_left"

                # Ey component
                p_mod_ey = figure(width=800, height=300, title=f"Modulus of Visibility (Ey) vs Time for Baseline {key}")
                p_mod_ey.line(time_points, moduli_over_time[key][:, 0, 1], line_width=2, legend_label=f"Ey Baseline {key}")
                p_mod_ey.xaxis.axis_label = 'Time (seconds)'
                p_mod_ey.yaxis.axis_label = 'Modulus of Visibility'
                p_mod_ey.legend.location = "top_left"

                # Phase components for Ex and Ey
                p_phase_ex = figure(width=800, height=300, title=f"Phase of Visibility (Ex) vs Time for Baseline {key}")
                p_phase_ex.line(time_points, phases_over_time[key][:, 0, 0], line_width=2, legend_label=f"Ex Baseline {key}")
                p_phase_ex.xaxis.axis_label = 'Time (seconds)'
                p_phase_ex.yaxis.axis_label = 'Phase of Visibility (radians)'
                p_phase_ex.legend.location = "top_left"

                p_phase_ey = figure(width=800, height=300, title=f"Phase of Visibility (Ey) vs Time for Baseline {key}")
                p_phase_ey.line(time_points, phases_over_time[key][:, 0, 1], line_width=2, legend_label=f"Ey Baseline {key}")
                p_phase_ey.xaxis.axis_label = 'Time (seconds)'
                p_phase_ey.yaxis.axis_label = 'Phase of Visibility (radians)'
                p_phase_ey.legend.location = "top_left"

                plots.extend([p_mod_ex, p_mod_ey, p_phase_ex, p_phase_ey])

            else:
                # No polarization, single component
                p_mod = figure(width=800, height=300, title=f"Modulus of Visibility vs Time for Baseline {key}")
                p_mod.line(time_points, moduli_over_time[key][:, 0], line_width=2, legend_label=f"Baseline {key}")
                p_mod.xaxis.axis_label = 'Time (seconds)'
                p_mod.yaxis.axis_label = 'Modulus of Visibility'
                p_mod.legend.location = "top_left"

                p_phase = figure(width=800, height=300, title=f"Phase of Visibility vs Time for Baseline {key}")
                p_phase.line(time_points, phases_over_time[key][:, 0], line_width=2, legend_label=f"Baseline {key}")
                p_phase.xaxis.axis_label = 'Time (seconds)'
                p_phase.yaxis.axis_label = 'Phase of Visibility (radians)'
                p_phase.legend.location = "top_left"

                plots.append(p_mod)
                plots.append(p_phase)

        return column(*plots, sizing_mode="stretch_both")

    else:
        # Matplotlib plotting
        num_baselines = len(baselines)
        fig1, ax1 = plt.subplots(num_baselines, 2, figsize=(15, 5 * num_baselines))
        
        for i, key in enumerate(baseline_keys):
            # If polarization is included, plot both Ex and Ey components
            if moduli_over_time[key].shape[-1] == 2:
                ax1[i, 0].plot(time_points, moduli_over_time[key][:, 0, 0], marker='o', label=f"Ex Baseline {key}")
                ax1[i, 0].plot(time_points, moduli_over_time[key][:, 0, 1], marker='o', label=f"Ey Baseline {key}")
                ax1[i, 0].set_xlabel('Time (seconds)')
                ax1[i, 0].set_ylabel('Modulus of Visibility')
                ax1[i, 0].set_title(f'Modulus of Visibility (Ex, Ey) vs Time for Baseline {key}')
                ax1[i, 0].legend()
                ax1[i, 0].grid(True)

                ax1[i, 1].plot(time_points, phases_over_time[key][:, 0, 0], marker='o', label=f"Ex Baseline {key}")
                ax1[i, 1].plot(time_points, phases_over_time[key][:, 0, 1], marker='o', label=f"Ey Baseline {key}")
                ax1[i, 1].set_xlabel('Time (seconds)')
                ax1[i, 1].set_ylabel('Phase of Visibility (radians)')
                ax1[i, 1].set_title(f'Phase of Visibility (Ex, Ey) vs Time for Baseline {key}')
                ax1[i, 1].legend()
                ax1[i, 1].grid(True)
            else:
                # No polarization, single component
                ax1[i, 0].plot(time_points, moduli_over_time[key][:, 0], marker='o', label=f"Baseline {key}")
                ax1[i, 0].set_xlabel('Time (seconds)')
                ax1[i, 0].set_ylabel('Modulus of Visibility')
                ax1[i, 0].set_title(f'Modulus of Visibility vs Time for Baseline {key}')
                ax1[i, 0].legend()
                ax1[i, 0].grid(True)

                ax1[i, 1].plot(time_points, phases_over_time[key][:, 0], marker='o', label=f"Baseline {key}")
                ax1[i, 1].set_xlabel('Time (seconds)')
                ax1[i, 1].set_ylabel('Phase of Visibility (radians)')
                ax1[i, 1].set_title(f'Phase of Visibility vs Time for Baseline {key}')
                ax1[i, 1].legend()
                ax1[i, 1].grid(True)

        plt.tight_layout()
        return fig1

def plot_heatmaps(moduli_over_time, phases_over_time, baselines, freqs, total_seconds, plotting='matplotlib'):
    """
    Create heatmaps for visibility modulus and phase using either Matplotlib or Bokeh.
    Returns appropriate plots based on the plotting library specified.
    Handles Ex and Ey polarization components if available.
    """
    baseline_keys = list(baselines.keys())

    if plotting == 'bokeh':
        plots = []
        for key in baseline_keys:
            if moduli_over_time[key].shape[-1] == 2:
                moduli_total = np.sqrt(moduli_over_time[key][:, :, 0]**2 + moduli_over_time[key][:, :, 1]**2)
                phases_total = np.sqrt(phases_over_time[key][:, :, 0]**2 + phases_over_time[key][:, :, 1]**2)
            else:
                moduli_total = moduli_over_time[key]
                phases_total = phases_over_time[key]

            p_mod = figure(width=800, height=300, title=f"Modulus of Visibility Heatmap for Baseline {key}")
            p_mod.image(image=[moduli_total], x=0, y=freqs[0] / 1e6, dw=total_seconds, dh=(freqs[-1] - freqs[0]) / 1e6, palette="Viridis256")
            p_mod.xaxis.axis_label = 'Time (seconds)'
            p_mod.yaxis.axis_label = 'Frequency (MHz)'

            p_phase = figure(width=800, height=300, title=f"Phase of Visibility Heatmap for Baseline {key}")
            p_phase.image(image=[phases_total], x=0, y=freqs[0] / 1e6, dw=total_seconds, dh=(freqs[-1] - freqs[0]) / 1e6, palette="Viridis256")
            p_phase.xaxis.axis_label = 'Time (seconds)'
            p_phase.yaxis.axis_label = 'Frequency (MHz)'

            plots.append(p_mod)
            plots.append(p_phase)

        return column(*plots, sizing_mode="stretch_both")

    else:
        num_baselines = len(baselines)
        fig2, ax2 = plt.subplots(num_baselines, 2, figsize=(15, 5 * num_baselines))

        for i, key in enumerate(baseline_keys):
            if moduli_over_time[key].shape[-1] == 2:
                moduli_total = np.sqrt(moduli_over_time[key][:, :, 0]**2 + moduli_over_time[key][:, :, 1]**2)
                phases_total = np.sqrt(phases_over_time[key][:, :, 0]**2 + phases_over_time[key][:, :, 1]**2)
            else:
                moduli_total = moduli_over_time[key]
                phases_total = phases_over_time[key]

            im0 = ax2[i, 0].imshow(moduli_total, aspect='auto', origin='lower',
                                   extent=[0, total_seconds, freqs[0] / 1e6, freqs[-1] / 1e6])
            ax2[i, 0].set_xlabel('Time (seconds)')
            ax2[i, 0].set_ylabel('Frequency (MHz)')
            ax2[i, 0].set_title(f'Modulus of Visibility for Baseline {key}')
            fig2.colorbar(im0, ax=ax2[i, 0])

            im1 = ax2[i, 1].imshow(phases_total, aspect='auto', origin='lower',
                                   extent=[0, total_seconds, freqs[0] / 1e6, freqs[-1] / 1e6])
            ax2[i, 1].set_xlabel('Time (seconds)')
            ax2[i, 1].set_ylabel('Frequency (MHz)')
            ax2[i, 1].set_title(f'Phase of Visibility for Baseline {key}')
            fig2.colorbar(im1, ax=ax2[i, 1])

        plt.tight_layout()
        return fig2

def plot_modulus_vs_frequency(moduli_over_time, baselines, freqs, time_index, plotting='matplotlib'):
    """
    Create the modulus of visibility vs frequency plot using either Matplotlib or Bokeh.
    Returns appropriate plots based on the plotting library specified.
    Handles Ex and Ey polarization components if available.
    """
    baseline_keys = list(baselines.keys())

    if plotting == 'bokeh':
        plots = []
        for key in baseline_keys:
            if moduli_over_time[key].shape[-1] == 2:
                moduli_total = np.sqrt(moduli_over_time[key][time_index, :, 0]**2 + moduli_over_time[key][time_index, :, 1]**2)
            else:
                moduli_total = moduli_over_time[key][time_index, :]

            p_mod = figure(width=800, height=300, title=f"Modulus of Visibility vs Frequency for Baseline {key} at Time {time_index} seconds")
            p_mod.line(freqs / 1e6, moduli_total, line_width=2, legend_label=f"Baseline {key}")
            p_mod.xaxis.axis_label = 'Frequency (MHz)'
            p_mod.yaxis.axis_label = 'Modulus of Visibility'
            p_mod.legend.location = "top_left"

            plots.append(p_mod)

        return column(*plots, sizing_mode="stretch_both")

    else:
        num_baselines = len(baselines)
        fig3, ax3 = plt.subplots(num_baselines, 1, figsize=(15, 5 * num_baselines))

        for i, key in enumerate(baseline_keys):
            if moduli_over_time[key].shape[-1] == 2:
                moduli_total = np.sqrt(moduli_over_time[key][time_index, :, 0]**2 + moduli_over_time[key][time_index, :, 1]**2)
            else:
                moduli_total = moduli_over_time[key][time_index, :]

            ax3[i].plot(freqs / 1e6, moduli_total, marker='o', label=f"Baseline {key}")
            ax3[i].set_xlabel('Frequency (MHz)')
            ax3[i].set_ylabel('Modulus of Visibility')
            ax3[i].set_title(f'Modulus of Visibility vs Frequency for Baseline {key} at Time {time_index} seconds')
            ax3[i].legend()
            ax3[i].grid(True)

        plt.tight_layout()
        return fig3
