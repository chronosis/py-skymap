from pathlib import Path

"""
Shared constants for skymap generation.

Centralising these here makes it easier to tweak behaviour without
digging through the main script.
"""

# Cache and output locations
CACHE_DB: Path = Path("gaia_cache/gaia_cache.db")
IMAGES_DIR: Path = Path("images")

# Gaia download defaults
CHUNK_SIZE: int = 1_000_000  # Stars per chunk (reduced to avoid timeouts)
DEFAULT_STAR_LIMIT: int = 50_000  # Default number of stars to download

# Gaia network / retry behaviour
# Previous defaults were:
# - max retries: 3 attempts
# - request timeout: 30 seconds
# - HEAD timeout: 10 seconds
# These values have been doubled to make the client more tolerant of
# transient Gaia service slowdowns.
GAIA_MAX_RETRIES: int = 6
GAIA_HTTP_TIMEOUT_SECONDS: int = 30
GAIA_HTTP_HEAD_TIMEOUT_SECONDS: int = 10
GAIA_MIN_MAGNITUDE_THRESHOLD: float = 16.5

# Distance / background classification
ASSUMED_BACKGROUND_DISTANCE_PC: float = 10_000_000.0  # 10 Mpc for effectively infinite distance
BRIGHT_BACKGROUND_ABS_MAG_THRESHOLD: float = 11.5
BACKGROUND_DISTANCE_THRESHOLD_PC: float = 100_000.0  # Stars beyond this are treated as background
MAX_DISTANCE_PC: float = 1e6  # General upper sanity bound for distances

# Visibility / plotting thresholds
VISIBLE_MAG_LIMIT: float = 11.5  # Faintest stars we render in the general field
VISIBLE_MAG_TRANSPARENCY_LIMIT: float = 7.0  # Stars brighter than this are fully opaque

# Point size scaling for general star field
STAR_POINT_MIN_SIZE: float = 3.0
STAR_POINT_MAX_SIZE: float = 132.0
STAR_POINT_MIN_ALPHA: float = 0.1
STAR_POINT_MAX_ALPHA: float = 0.9

# Figure / image dimensions (used for proportional text sizing)
FIGURE_SIZE_INCHES: float = 24.0
FIGURE_DPI: int = 150

# Label font sizes: proportional to figure size (fontsize = ratio * FIGURE_SIZE_INCHES).
# Bumped from original 7/8/9 to improve readability.
LABEL_FONT_SIZE_RATIO_SMALL: float = 12.0 / 24.0   # Bright star labels
LABEL_FONT_SIZE_RATIO_MED: float = 14.0 / 24.0     # DSO labels, axis tick labels
LABEL_FONT_SIZE_RATIO_LARGE: float = 16.0 / 24.0  # Galaxy, Magellanic cloud, Sgr A*, Sol

# Text stroke: semi-transparent black outline for legibility on varied backgrounds
TEXT_STROKE_LINEWIDTH: float = 4.0 * 72.0 / 150.0  # ~2 pixels at FIGURE_DPI
TEXT_STROKE_COLOR: str = "black"
TEXT_STROKE_ALPHA: float = 0.9

