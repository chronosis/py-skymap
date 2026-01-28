import numpy as np
import matplotlib.pyplot as plt

from .math3d import transform_to_target_frame, get_relative_coords


def bp_rp_to_rgb(bp_rp, alpha=0.3):
    """Convert Gaia BP-RP color index to RGB color with optional transparency."""
    bp_rp = np.asarray(bp_rp)

    # Handle NaN values (stars without color data) - return white
    nan_mask = ~np.isfinite(bp_rp)

    # Clamp BP-RP to reasonable range for color mapping
    bp_rp_clamped = np.clip(bp_rp, -0.5, 3.0)

    # Normalize to [0, 1] for colormap
    bp_rp_normalized = (bp_rp_clamped + 0.5) / 3.5

    colors = np.zeros((len(bp_rp) if bp_rp.ndim > 0 else 1, 4))

    # Blue stars (bp_rp < 0.5): blue to cyan
    blue_mask = bp_rp_normalized < 0.2
    if np.any(blue_mask):
        colors[blue_mask, 0] = 0.0  # R
        colors[blue_mask, 1] = 0.5 + (bp_rp_normalized[blue_mask] / 0.2) * 0.5
        colors[blue_mask, 2] = 1.0  # B

    # White/neutral stars (0.5 < bp_rp < 1.5): white to yellow
    neutral_mask = (bp_rp_normalized >= 0.2) & (bp_rp_normalized < 0.5)
    if np.any(neutral_mask):
        neutral_norm = (bp_rp_normalized[neutral_mask] - 0.2) / (0.5 - 0.2)
        colors[neutral_mask, 0] = 0.5 + neutral_norm * 0.5
        colors[neutral_mask, 1] = 1.0
        colors[neutral_mask, 2] = 1.0 - neutral_norm * 0.5

    # Orange/red stars (bp_rp > 1.5): yellow to red
    red_mask = bp_rp_normalized >= 0.5
    if np.any(red_mask):
        colors[red_mask, 0] = 1.0
        red_norm = (bp_rp_normalized[red_mask] - 0.5) / (1.0 - 0.5)
        colors[red_mask, 1] = 1.0 - red_norm
        colors[red_mask, 2] = 0.0

    colors[:, 3] = alpha
    colors[nan_mask] = [1.0, 1.0, 1.0, alpha]

    # Safety clamp
    colors = np.clip(colors * 2, 0.0, 1.0)

    if np.isscalar(bp_rp):
        return tuple(colors[0])
    return colors


def calculate_point_size_by_magnitude(magnitude, min_size=1, max_size=10, magnitude_brightest=10):
    """Calculate point size based on apparent magnitude using linear interpolation."""
    magnitude = np.asarray(magnitude)

    point_sizes = np.zeros_like(magnitude, dtype=float)

    above_brightest = magnitude > magnitude_brightest
    point_sizes[above_brightest] = 0.0

    below_zero = magnitude < 0
    point_sizes[below_zero] = max_size

    in_range = (magnitude >= 0) & (magnitude <= magnitude_brightest)
    if np.any(in_range):
        mag_range = magnitude_brightest - 0.0
        size_range = max_size - min_size
        point_sizes[in_range] = max_size - magnitude[in_range] * size_range / mag_range

    if np.isscalar(magnitude):
        return float(point_sizes)
    return point_sizes


# The remaining plotting helpers (bright objects, galaxies, DSOs, Sol, etc.)
# are left in skymap-gen.py for now but can be progressively migrated here
# by moving their function bodies and importing them from this module.

