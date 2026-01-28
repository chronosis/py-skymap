import numpy as np
import matplotlib.pyplot as plt

from .math3d import transform_to_target_frame, get_relative_coords


def bp_rp_to_rgb(bp_rp, alpha=0.3):
    """Convert Gaia BP-RP color index to RGB(A).

    `alpha` can be a scalar (uniform transparency) or an array matching
    the shape of `bp_rp` (per-star transparency).
    """
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

    # Handle alpha as scalar or array
    alpha_arr = np.asarray(alpha)
    if alpha_arr.ndim == 0:
        colors[:, 3] = float(alpha_arr)
        nan_alpha = float(alpha_arr)
    else:
        # Assume shape is broadcast-compatible with bp_rp
        colors[:, 3] = alpha_arr
        nan_alpha = np.broadcast_to(alpha_arr, colors.shape[0])[0]

    # NaN colors: white with same alpha
    colors[nan_mask] = [1.0, 1.0, 1.0, nan_alpha]

    # Make star-field colours softer and less saturated by blending toward white.
    colors[:, :3] = np.clip(colors[:, :3], 0.0, 1.0)
    colors[:, :3] = 0.5 * colors[:, :3] + 0.5  # move halfway toward white
    colors[:, 3] = np.clip(colors[:, 3], 0.0, 1.0)

    if np.isscalar(bp_rp):
        return tuple(colors[0])
    return colors


def calculate_point_size_by_magnitude(magnitude, min_size=1, max_size=10, magnitude_brightest=10):
    """Calculate point size based on apparent magnitude using piecewise interpolation.
    
    Uses a piecewise function: bright stars (mag <= 7.0) get the full dynamic range,
    while faint stars (mag > 7.0) are clamped to min_size. This ensures distant
    background stars stay small even when max_size is increased.
    """
    magnitude = np.asarray(magnitude)

    point_sizes = np.zeros_like(magnitude, dtype=float)

    above_brightest = magnitude > magnitude_brightest
    point_sizes[above_brightest] = 0.0

    below_zero = magnitude < 0
    point_sizes[below_zero] = max_size

    mag_range = magnitude_brightest - 0.0
    size_range = max_size - min_size
    in_range = (magnitude >= 0) & (magnitude <= magnitude_brightest)
    if np.any(in_range):
        # Use a piecewise function: bright stars get dynamic range, faint stars stay near min_size
        # - Bright stars (mag <= 7.0): use full dynamic range from max_size down
        # - Faint stars (mag > 7.0): clamp to min_size (or very close) so they don't scale with max_size
        mag_threshold = 7.0
        
        # Create masks directly on the full magnitude array
        bright_mask = in_range & (magnitude <= mag_threshold)
        faint_mask = in_range & (magnitude > mag_threshold)
        
        if np.any(bright_mask):
            # Bright stars: linear mapping from mag 0-7 to size max_size down to min_size
            bright_mag_norm = magnitude[bright_mask] / mag_threshold  # 0 to 1
            point_sizes[bright_mask] = max_size - bright_mag_norm * size_range
        
        if np.any(faint_mask):
            # Faint stars: all get min_size (or very close) regardless of max_size
            # This ensures distant background stars stay small even when max_size increases
            point_sizes[faint_mask] = min_size
    
    if np.isscalar(magnitude):
        return float(point_sizes)
    return point_sizes


# The remaining plotting helpers (bright objects, galaxies, DSOs, Sol, etc.)
# are left in skymap-gen.py for now but can be progressively migrated here
# by moving their function bodies and importing them from this module.

