# src/baseline.py
import sys
import numpy as np

def generate_baselines(antennas):
    """
    Generate baselines based on antenna positions, using Numbers as keys.

    Parameters:
    antennas (dict): Dictionary of antenna metadata and positions.

    Returns:
    dict: Dictionary of baselines between antennas with Numbers as keys.
    """
    baselines = {}
    # Extract antenna numbers and positions
    antenna_numbers = {ant["Number"]: np.array(ant["Position"]) for ant in antennas.values()}

    # Sort antenna numbers to ensure all combinations are considered
    sorted_antenna_numbers = sorted(antenna_numbers.keys())

    # Generate baselines
    for i, ant1 in enumerate(sorted_antenna_numbers):
        for ant2 in sorted_antenna_numbers[i:]:  # Ensure ant1 <= ant2
            pos1 = antenna_numbers[ant1]
            pos2 = antenna_numbers[ant2]
            baselines[(ant1, ant2)] = pos2 - pos1
    # Calculate total memory usage in MB
    total_memory_bytes = sys.getsizeof(baselines) + sum(
        sys.getsizeof(key) + sys.getsizeof(value) for key, value in baselines.items()
    )
    total_memory_mb = total_memory_bytes / (1024 * 1024)
    print(f"Total memory used by baselines: {total_memory_mb:.4f} MB")

    return baselines

