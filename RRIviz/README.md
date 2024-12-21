
# Radio Astronomy Visibility Simulator

## Overview

The **Radio Astronomy Visibility Simulator** is a Python-based application designed for simulating and visualizing complex visibilities in radio astronomy. The simulator calculates visibility data using antenna positions, baselines, sky models, and source catalogs like GLEAM. It supports multiple plotting libraries and provides options for polarized and non-polarized visibility calculations.

## Features

- **Antenna Configuration**: Generate antenna positions or load them from a file.
- **Baseline Calculation**: Dynamically compute baselines between antennas.
- **Sky Model Integration**: Utilize global and GLEAM sky models for accurate source data.
- **Visibility Calculation**:
  - Optimized and original methods for visibility computation.
  - Support for polarized visibility.
- **Interactive Plots**: Visualize results using Matplotlib or Bokeh.
- **Customizability**: Adjust parameters like field of view, frequency, and antenna spacing.
- **Data Persistence**: Save visibility data and skymodel plots for further analysis.

---

## Installation

1. Clone the repository:
   ```bash
   git clone <repository_url>
   cd <repository_directory>
   ```

2. Create a virtual environment and activate it:
   ```bash
   python -m venv matvis-env
   source matvis-env/bin/activate  # On Windows, use matvis-env\Scripts\activate
   ```

3. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```

---

## Usage

### Command-Line Options
Run the simulator with various configurations using the command line:

```bash
python src/main.py [OPTIONS]
```

Key arguments:
- `--num_antennas`: Number of antennas to simulate (default: 3).
- `--spacing`: Spacing between antennas in meters (default: 14.0).
- `--lat`, `--lon`, `--height`: Latitude, longitude, and height of the observation location.
- `--antenna_positions_file`: File path for custom antenna positions.
- `--use_gleam`: Use the GLEAM catalog for source modeling.
- `--gleam_flux_limit`: Minimum flux of sources in Jy (default: 1.0).
- `--plotting`: Choose between `matplotlib` and `bokeh` for visualization.
- `--use_polarization`: Enable polarized visibility calculation.
- `--plot_skymodel_every_hour`: Generate and plot sky model hourly.

### Example
Simulate visibility with 5 antennas, using GLEAM sources and Bokeh for plotting:
```bash
python src/main.py --num_antennas 5 --use_gleam --gleam_flux_limit 2.0 --plotting bokeh
```

---

## File Structure

```
.
├── README.md                # Project documentation
├── data
│   └── antenna.txt          # Sample antenna positions file
├── requirements.txt         # Python dependencies
├── setup.py                 # Package setup
├── src                      # Source code
│   ├── antenna.py           # Antenna position handling
│   ├── baseline.py          # Baseline generation
│   ├── main.py              # Entry point for the simulator
│   ├── observation.py       # Observation location and time
│   ├── plot.py              # Plotting functions
│   ├── skymodel.py          # Sky model integration
│   ├── source.py            # Source modeling (GLEAM and test sources)
│   ├── visibility.py        # Visibility calculations
│   └── visibility_results.json # Output data
└── tests                    # Unit tests
    ├── test_antenna.py
    ├── test_baseline.py
    ├── test_main.py
    └── test_visibility.py
```

---

## Development

### Requirements
- Python 3.8+
- Key libraries: `numpy`, `astropy`, `healpy`, `matplotlib`, `bokeh`, `astroquery`.

### Running Tests
Run unit tests using `pytest`:
```bash
pytest tests/
```

---

## Contributors

- **Kartik Mandar** - Developer
- Open-source contributions are welcome! Fork the repo and submit a pull request.

---

## License

This project is licensed under the MIT License. See `LICENSE` for details.
