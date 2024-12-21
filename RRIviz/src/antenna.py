# src/antenna.py
import sys

def read_antenna_positions(file_path):
    """
    Reads antenna positions and metadata from a file.

    The file should contain antenna information in the format:
    Name  Number  BeamID  E  N  U

    Parameters:
    file_path (str): Path to the antenna position file.

    Returns:
    dict: Dictionary with antenna indices as keys and dictionaries of metadata and positions as values.

    Raises:
    ValueError: If the file is empty or has invalid data.
    """
    antennas = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Ensure the file has more than just the header
    if len(lines) <= 1:
        raise ValueError("The antenna positions file is empty or has no valid data.")

    for i, line in enumerate(lines):
        # Skip the header line or empty lines
        if i == 0 or not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 6:
            raise ValueError(f"Invalid antenna position in line {i+1}: {line}")
        # Extract metadata and positions
        try:
            name = parts[0]
            number = int(parts[1])
            beam_id = int(parts[2])
            e, n, u = map(float, parts[3:6])
            antennas[i - 1] = {
                "Name": name,
                "Number": number,
                "BeamID": beam_id,
                "Position": (e, n, u)
            }
        except ValueError:
            raise ValueError(f"Could not parse data in line {i+1}: {line}")
        
    # Calculate total memory usage in MB
    total_memory_bytes = sys.getsizeof(antennas) + sum(
        sys.getsizeof(value) + sum(sys.getsizeof(v) for v in value.values()) 
        for value in antennas.values()
    )
    total_memory_mb = total_memory_bytes / (1024 * 1024)
    print(f"Total memory used by antennas: {total_memory_mb:.4f} MB")

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
        antennas[i] = (i * spacing, 0, 0)  # Antennas aligned along the EW axis

    return antennas
