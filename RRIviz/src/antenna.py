# src/antenna.py

def read_antenna_positions(file_path):
    """
    Reads antenna positions from a file.

    The file should contain one antenna position per line in the format:
    x y z

    Returns:
    dict: Dictionary with antenna indices and their positions.
    """
    antennas = {}
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) != 3:
                raise ValueError(f"Invalid antenna position in line {i+1}: {line}")
            x, y, z = map(float, parts)
            antennas[i] = (x, y, z)
    return antennas

def generate_antennas(num_antennas=3, spacing=14.0):
    """
    Generate antenna positions in an east-west line.

    Parameters:
    num_antennas (int): Number of antennas to generate.
    spacing (float): Distance between each antenna in meters.

    Returns:
    dict: Dictionary with antenna indices and their positions.
    """
    antennas = {}

    for i in range(num_antennas):
        antennas[i] = (i * spacing, 0, 0)  # Antennas aligned along the x-axis

    return antennas
