from bokeh.plotting import figure, show
from bokeh.layouts import column
import matplotlib.pyplot as plt

def plot_visibility(moduli_over_time, phases_over_time, baselines, time_points, freqs, total_seconds, plotting='matplotlib'):
    """
    Create the modulus and phase of visibility vs time using either Matplotlib or Bokeh.
    Returns appropriate plots based on the plotting library specified.
    """
    baseline_keys = list(baselines.keys())

    if plotting == 'bokeh':
        # Bokeh plotting
        plots = []
        for key in baseline_keys:
            # Create plot for Modulus of Visibility vs Time
            p_mod = figure(width=800, height=300, title=f"Modulus of Visibility vs Time for Baseline {key}")
            p_mod.line(time_points, moduli_over_time[key][:, 0], line_width=2, legend_label=f"Baseline {key}")
            p_mod.xaxis.axis_label = 'Time (seconds)'
            p_mod.yaxis.axis_label = 'Modulus of Visibility'
            p_mod.legend.location = "top_left"

            # Create plot for Phase of Visibility vs Time
            p_phase = figure(width=800, height=300, title=f"Phase of Visibility vs Time for Baseline {key}")
            p_phase.line(time_points, phases_over_time[key][:, 0], line_width=2, legend_label=f"Baseline {key}")
            p_phase.xaxis.axis_label = 'Time (seconds)'
            p_phase.yaxis.axis_label = 'Phase of Visibility (radians)'
            p_phase.legend.location = "top_left"

            plots.append(p_mod)
            plots.append(p_phase)

        # Return Bokeh scrollable layout
        return column(*plots, sizing_mode="stretch_both")

    else:
        # Matplotlib plotting
        num_baselines = len(baselines)
        fig1, ax1 = plt.subplots(num_baselines, 2, figsize=(15, 5 * num_baselines))
        
        for i, key in enumerate(baseline_keys):
            # Plot modulus of visibility vs time
            ax1[i, 0].plot(time_points, moduli_over_time[key][:, 0], marker='o', label=f"Baseline {key}")
            ax1[i, 0].set_xlabel('Time (seconds)')
            ax1[i, 0].set_ylabel('Modulus of Visibility')
            ax1[i, 0].set_title(f'Modulus of Visibility vs Time for Baseline {key}')
            ax1[i, 0].legend()
            ax1[i, 0].grid(True)

            # Plot phase of visibility vs time
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
    """
    baseline_keys = list(baselines.keys())

    if plotting == 'bokeh':
        plots = []
        for key in baseline_keys:
            # Create heatmap for Modulus of Visibility
            p_mod = figure(width=800, height=300, title=f"Modulus of Visibility Heatmap for Baseline {key}")
            p_mod.image(image=[moduli_over_time[key]], x=0, y=freqs[0] / 1e6, dw=total_seconds, dh=(freqs[-1] - freqs[0]) / 1e6, palette="Viridis256")
            p_mod.xaxis.axis_label = 'Time (seconds)'
            p_mod.yaxis.axis_label = 'Frequency (MHz)'

            # Create heatmap for Phase of Visibility
            p_phase = figure(width=800, height=300, title=f"Phase of Visibility Heatmap for Baseline {key}")
            p_phase.image(image=[phases_over_time[key]], x=0, y=freqs[0] / 1e6, dw=total_seconds, dh=(freqs[-1] - freqs[0]) / 1e6, palette="Viridis256")
            p_phase.xaxis.axis_label = 'Time (seconds)'
            p_phase.yaxis.axis_label = 'Frequency (MHz)'

            plots.append(p_mod)
            plots.append(p_phase)

        return column(*plots, sizing_mode="stretch_both")

    else:
        # Matplotlib plotting
        num_baselines = len(baselines)
        fig2, ax2 = plt.subplots(num_baselines, 2, figsize=(15, 5 * num_baselines))

        for i, key in enumerate(baseline_keys):
            # Heatmap for modulus of visibility
            im0 = ax2[i, 0].imshow(moduli_over_time[key], aspect='auto', origin='lower',
                                   extent=[0, total_seconds, freqs[0] / 1e6, freqs[-1] / 1e6])
            ax2[i, 0].set_xlabel('Time (seconds)')
            ax2[i, 0].set_ylabel('Frequency (MHz)')
            ax2[i, 0].set_title(f'Modulus of Visibility for Baseline {key}')
            fig2.colorbar(im0, ax=ax2[i, 0])

            # Heatmap for phase of visibility
            im1 = ax2[i, 1].imshow(phases_over_time[key], aspect='auto', origin='lower',
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
    """
    baseline_keys = list(baselines.keys())

    if plotting == 'bokeh':
        plots = []
        for key in baseline_keys:
            # Create plot for Modulus of Visibility vs Frequency
            p_mod = figure(width=800, height=300, title=f"Modulus of Visibility vs Frequency for Baseline {key} at Time {time_index} seconds")
            p_mod.line(freqs / 1e6, moduli_over_time[key][time_index, :], line_width=2, legend_label=f"Baseline {key}")
            p_mod.xaxis.axis_label = 'Frequency (MHz)'
            p_mod.yaxis.axis_label = 'Modulus of Visibility'
            p_mod.legend.location = "top_left"

            plots.append(p_mod)

        return column(*plots, sizing_mode="stretch_both")

    else:
        # Matplotlib plotting
        num_baselines = len(baselines)
        fig3, ax3 = plt.subplots(num_baselines, 1, figsize=(15, 5 * num_baselines))

        for i, key in enumerate(baseline_keys):
            # Plot modulus of visibility vs frequency for a specific time index
            ax3[i].plot(freqs / 1e6, moduli_over_time[key][time_index, :], marker='o', label=f"Baseline {key}")
            ax3[i].set_xlabel('Frequency (MHz)')
            ax3[i].set_ylabel('Modulus of Visibility')
            ax3[i].set_title(f'Modulus of Visibility vs Frequency for Baseline {key} at Time {time_index} seconds')
            ax3[i].legend()
            ax3[i].grid(True)

        plt.tight_layout()
        return fig3
