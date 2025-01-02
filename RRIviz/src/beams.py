# beams.py

import numpy as np

def calculate_A_theta(theta, theta_HPBW):
    """
    Calculate the Gaussian primary beam pattern A(theta).

    Parameters:
    theta (ndarray): Zenith angles in radians.
    theta_HPBW (float): Half Power Beam Width (HPBW) in radians.

    Returns:
    ndarray: The primary beam pattern A(theta) evaluated at theta.
    """
    return np.exp(- (theta / (np.sqrt(2) * theta_HPBW))**2 ) ** 2
