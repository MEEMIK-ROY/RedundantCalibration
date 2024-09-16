# src/baseline.py

import numpy as np

def generate_baselines(antennas):
    """
    Generate baselines based on antenna positions.

    Parameters:
    antennas (dict): Dictionary of antenna positions.

    Returns:
    dict: Dictionary of baselines between antennas.
    """
    baselines = {}
    antenna_indices = list(antennas.keys())

    for i in range(len(antenna_indices)):
        for j in range(i + 1, len(antenna_indices)):
            ant1 = antenna_indices[i]
            ant2 = antenna_indices[j]
            pos1 = np.array(antennas[ant1])
            pos2 = np.array(antennas[ant2])
            baseline = pos2 - pos1
            baselines[(ant1, ant2)] = baseline

    return baselines
