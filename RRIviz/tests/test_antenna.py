# tests/test_antenna.py

from src.antenna import antennas

def test_antenna_positions():
    # Test the number of antennas
    assert len(antennas) == 3
    # Test individual antenna positions
    assert antennas[0] == (0, 0, 0)
    assert antennas[1] == (14, 0, 0)
    assert antennas[2] == (28, 0, 0)
