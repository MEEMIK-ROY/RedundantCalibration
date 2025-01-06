from astropy.io import fits
import os

# Define the relative file path
relative_path = '../data/NF_HERA_Vivaldi_efield_beam_healpix.fits'

# Compute the absolute path
file_path = os.path.abspath(relative_path)
print("Attempting to open file at:", file_path)

# Try to open the file
try:
    with fits.open(file_path) as hdul:
        # Print a summary of the file
        hdul.info()
        
        # Access specific HDU data (e.g., first HDU)
        data = hdul[0].data
        header = hdul[0].header

        # Print the header
        print("\nHeader:")
        print(repr(header))

        # If data exists, print its shape
        if data is not None:
            print("\nData shape:", data.shape)
            print("Data sample:", data[0])
except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
