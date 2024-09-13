# src/plot.py

import matplotlib.pyplot as plt

def plot_visibility(moduli_over_time, phases_over_time, baselines, time_points, freqs, total_seconds):
    """
    Create the modulus and phase of visibility vs time and return the figures for combined display.
    """
    fig1, ax1 = plt.subplots(3, 2, figsize=(15, 15))
    baseline_keys = list(baselines.keys())

    for i, key in enumerate(baseline_keys):
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

    return fig1

def plot_heatmaps(moduli_over_time, phases_over_time, baselines, freqs, total_seconds):
    """
    Create the heatmaps for visibility and return the figure.
    """
    fig2, ax2 = plt.subplots(3, 2, figsize=(15, 15))
    baseline_keys = list(baselines.keys())

    for i, key in enumerate(baseline_keys):
        im0 = ax2[i, 0].imshow(moduli_over_time[key], aspect='auto', origin='lower', extent=[0, total_seconds, freqs[0] / 1e6, freqs[-1] / 1e6])
        ax2[i, 0].set_xlabel('Time (seconds)')
        ax2[i, 0].set_ylabel('Frequency (MHz)')
        ax2[i, 0].set_title(f'Modulus of Visibility for Baseline {key}')
        fig2.colorbar(im0, ax=ax2[i, 0])

        im1 = ax2[i, 1].imshow(phases_over_time[key], aspect='auto', origin='lower', extent=[0, total_seconds, freqs[0] / 1e6, freqs[-1] / 1e6])
        ax2[i, 1].set_xlabel('Time (seconds)')
        ax2[i, 1].set_ylabel('Frequency (MHz)')
        ax2[i, 1].set_title(f'Phase of Visibility for Baseline {key}')
        fig2.colorbar(im1, ax=ax2[i, 1])

    return fig2

def plot_modulus_vs_frequency(moduli_over_time, baselines, freqs, time_index):
    """
    Create the modulus of visibility vs frequency plot and return the figure.
    """
    fig3, ax3 = plt.subplots(3, 1, figsize=(15, 15))
    baseline_keys = list(baselines.keys())

    for i, key in enumerate(baseline_keys):
        ax3[i].plot(freqs / 1e6, moduli_over_time[key][time_index, :], marker='o', label=f"Baseline {key}")
        ax3[i].set_xlabel('Frequency (MHz)')
        ax3[i].set_ylabel('Modulus of Visibility')
        ax3[i].set_title(f'Modulus of Visibility vs Frequency for Baseline {key} at Time {time_index} seconds')
        ax3[i].legend()
        ax3[i].grid(True)

    return fig3
