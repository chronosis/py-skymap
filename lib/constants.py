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
CHUNK_SIZE: int = 500_000  # Stars per chunk (reduced to avoid timeouts)
DEFAULT_STAR_LIMIT: int = 50_000  # Default number of stars to download

# Distance / background classification
ASSUMED_BACKGROUND_DISTANCE_PC: float = 10_000_000.0  # 10 Mpc for effectively infinite distance
BRIGHT_BACKGROUND_ABS_MAG_THRESHOLD: float = 10.0
BACKGROUND_DISTANCE_THRESHOLD_PC: float = 100_000.0  # Stars beyond this are treated as background
MAX_DISTANCE_PC: float = 1e6  # General upper sanity bound for distances

# Visibility / plotting thresholds
VISIBLE_MAG_LIMIT: float = 10.5  # Faintest stars we render in the general field

# Point size scaling for general star field
STAR_POINT_MIN_SIZE: float = 1.5
STAR_POINT_MAX_SIZE: float = 16.0

