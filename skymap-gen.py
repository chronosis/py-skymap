import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.patheffects as patheffects
import os
import sys
import time
import argparse
import sqlite3
from pathlib import Path
from astropy.coordinates import SkyCoord, Distance, CartesianRepresentation
from astropy import units as u
from astroquery.gaia import Gaia
from astropy.table import Table

from lib.progress import HAS_TQDM, tqdm

from lib.constants import (
    CACHE_DB,
    IMAGES_DIR,
    CHUNK_SIZE,
    DEFAULT_STAR_LIMIT,
    ASSUMED_BACKGROUND_DISTANCE_PC,
    BRIGHT_BACKGROUND_ABS_MAG_THRESHOLD,
    BACKGROUND_DISTANCE_THRESHOLD_PC,
    MAX_DISTANCE_PC,
    VISIBLE_MAG_LIMIT,
    VISIBLE_MAG_TRANSPARENCY_LIMIT,
    STAR_POINT_MIN_SIZE,
    STAR_POINT_MAX_SIZE,
    STAR_POINT_MIN_ALPHA,
    STAR_POINT_MAX_ALPHA,
    FIGURE_SIZE_INCHES,
    LABEL_FONT_SIZE_RATIO_SMALL,
    LABEL_FONT_SIZE_RATIO_MED,
    LABEL_FONT_SIZE_RATIO_LARGE,
    TEXT_STROKE_LINEWIDTH,
    TEXT_STROKE_COLOR,
    TEXT_STROKE_ALPHA,
)
from lib.sqlite_helper import (
    init_database as _sqlite_init_database,
    get_star_count as _sqlite_get_star_count,
    get_target_star_from_cache as _sqlite_get_target_star_from_cache,
    cache_target_star as _sqlite_cache_target_star,
    clear_star_positions_for_target as _sqlite_clear_star_positions_for_target,
    insert_star_positions_batch as _sqlite_insert_star_positions_batch,
    load_stars_from_cache as _sqlite_load_stars_from_cache,
    check_sky_coverage_bias as _sqlite_check_sky_coverage_bias,
)
from lib.math3d import (
    get_relative_coords as _math_get_relative_coords,
    transform_to_target_frame as _math_transform_to_target_frame,
)
from lib.plot_helper import (
    bp_rp_to_rgb as _plot_bp_rp_to_rgb,
    calculate_point_size_by_magnitude as _plot_calculate_point_size_by_magnitude,
)
from lib.star_data import (
    get_bright_galaxies as _stardata_get_bright_galaxies,
    get_bright_stars as _stardata_get_bright_stars,
    get_bright_deep_sky_objects as _stardata_get_bright_deep_sky_objects,
    get_magellanic_clouds as _stardata_get_magellanic_clouds,
    get_known_target_parallax_override as _stardata_get_known_target_parallax_override,
)


def _label_font_size(ratio):
    """Font size proportional to figure dimensions."""
    return ratio * FIGURE_SIZE_INCHES


def _text_stroke_effects():
    """Semi-transparent black 1-pixel stroke for text legibility."""
    stroke_color = (*mcolors.to_rgb(TEXT_STROKE_COLOR), TEXT_STROKE_ALPHA)
    return [
        patheffects.withStroke(linewidth=TEXT_STROKE_LINEWIDTH, foreground=stroke_color),
        patheffects.Normal(),
    ]


def bp_rp_to_rgb(bp_rp, alpha=0.3):
    """Delegate to shared color helper."""
    return _plot_bp_rp_to_rgb(bp_rp, alpha=alpha)

def init_database(db_path):
    """Initialize or return SQLite database connection via shared helper."""
    return _sqlite_init_database(db_path)

def get_star_count(db_path):
    """Get the number of stars in the database via shared helper."""
    return _sqlite_get_star_count(db_path)


def get_target_star_from_cache(db_path, target_star_name):
    """Check if target star is in cache via shared helper."""
    return _sqlite_get_target_star_from_cache(db_path, target_star_name)

def get_target_star_from_gaia(target_star_name):
    """Query Gaia database for target star by name using crossmatch.
    
    Returns a dict with source_id, ra, dec, parallax, phot_g_mean_mag,
    or None if not found.
    """
    print(f"  Querying Gaia for target star: {target_star_name}...")
    
    try:
        # First, get coordinates from Simbad
        simbad_coord = SkyCoord.from_name(target_star_name)
        ra_deg = simbad_coord.ra.deg
        dec_deg = simbad_coord.dec.deg
        
        print(f"  Found coordinates from Simbad: RA={ra_deg:.6f}°, Dec={dec_deg:.6f}°")
        
        # Query Gaia for stars near this position with a small search radius
        # Use a 1 arcminute search radius to find the target star
        search_radius_arcsec = 60  # 1 arcminute
        query = f"""
        SELECT TOP 10
            source_id, ra, dec, parallax, phot_g_mean_mag
        FROM gaiadr3.gaia_source
        WHERE 1=CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra_deg}, {dec_deg}, {search_radius_arcsec / 3600.0})
        )
        ORDER BY phot_g_mean_mag ASC
        """
        
        job = Gaia.launch_job(query)
        results = job.get_results()
        
        if len(results) == 0:
            print(f"  Warning: No stars found in Gaia near {target_star_name}")
            return None
        
        # Use the brightest star (lowest magnitude) as the target
        target = results[0]
        print(f"  Found target star in Gaia:")
        print(f"    Source ID: {target['source_id']}")
        print(f"    RA: {target['ra']:.6f}°, Dec: {target['dec']:.6f}°")
        print(f"    Parallax: {target['parallax']:.4f} mas")
        print(f"    Magnitude: {target['phot_g_mean_mag']:.2f}")
        
        return {
            'source_id': int(target['source_id']),
            'ra': float(target['ra']),
            'dec': float(target['dec']),
            'parallax': float(target['parallax']),
            'phot_g_mean_mag': float(target['phot_g_mean_mag'])
        }
        
    except Exception as e:
        print(f"  Error querying Gaia for target star: {e}")
        return None

def cache_target_star(db_path, target_data):
    """Cache the target star via shared SQLite helper."""
    return _sqlite_cache_target_star(db_path, target_data)


def clear_star_positions_for_target(db_path, target_star_name):
    """Remove all star_positions_3d rows for the given target via shared helper."""
    return _sqlite_clear_star_positions_for_target(db_path, target_star_name)


def insert_star_positions_batch(db_path, target_star_name, rows):
    """Insert a batch of positions via shared SQLite helper."""
    return _sqlite_insert_star_positions_batch(db_path, target_star_name, rows)


def load_stars_from_cache(db_path, limit=None, offset=0):
    """Load stars from SQLite cache via shared helper (no Gaia queries)."""
    return _sqlite_load_stars_from_cache(db_path, limit=limit, offset=offset)

def check_sky_coverage_bias(db_path):
    """Check if database has biased sky coverage (e.g., only northern hemisphere).
    Returns (has_bias, message) where has_bias is True if bias detected."""
    if not db_path.exists():
        return False, None
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    # Check dec distribution
    cursor.execute("SELECT COUNT(*) FROM gaia_source WHERE dec < 0")
    n_neg = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM gaia_source")
    total = cursor.fetchone()[0]
    conn.close()
    
    if total == 0:
        return False, None
    
    neg_ratio = n_neg / total
    # If less than 10% have dec < 0, likely biased (should be ~50% for uniform sky)
    if neg_ratio < 0.1:
        return True, (
            f"WARNING: Database appears to have biased sky coverage. "
            f"Only {n_neg:,} stars ({100*neg_ratio:.1f}%) have dec < 0 (should be ~50%). "
            f"This is likely due to ORDER BY source_id in old downloads. "
            f"Re-download with --force-refresh to get uniform coverage."
        )
    
    return False, None

def ensure_cache_populated(force_refresh=False, star_limit=None):
    """Ensure the SQLite cache is populated with Gaia data.
    
    This function checks the cache first and only queries Gaia API if:
    - Cache is empty, OR
    - force_refresh=True, OR
    - Cache has fewer stars than star_limit
    
    All processing operations should use load_stars_from_cache() to read from
    the SQLite database, never querying Gaia directly.
    
    Args:
        force_refresh: If True, delete existing cache and re-download
        star_limit: Maximum number of stars to cache (None = no limit)
    
    Returns:
        Number of stars in the cache
    """
    # Initialize database
    if force_refresh and CACHE_DB.exists():
        print(f"Removing existing database {CACHE_DB}...")
        CACHE_DB.unlink()
    
    # Check if cache already exists and is sufficient
    existing_count = get_star_count(CACHE_DB)
    if existing_count > 0 and not force_refresh:
        # Check for sky coverage bias
        has_bias, bias_msg = check_sky_coverage_bias(CACHE_DB)
        if has_bias:
            print(f"\n{bias_msg}\n")
        if star_limit is not None and existing_count >= star_limit:
            print(f"Cache has {existing_count:,} stars (meets limit of {star_limit:,}). Using cached data.")
            return existing_count
        else:
            print(f"Cache has {existing_count:,} stars. Checking if more needed...")
            if star_limit is None:
                # No limit specified, cache is sufficient
                print("Using existing cache (no limit specified).")
                return existing_count
            else:
                # Need more stars, continue downloading from Gaia then VizieR
                print(f"Continuing download to reach {star_limit:,} stars...")
                _download_from_gaia(star_limit, existing_count)
                _download_from_vizier(star_limit)
                return get_star_count(CACHE_DB)
    
    # Cache is empty or force_refresh, download from Gaia then VizieR
    _download_from_gaia(star_limit, 0)
    _download_from_vizier(star_limit)
    return get_star_count(CACHE_DB)

def _update_progress_bar(pbar, chunk_num, chunk_inserted, duplicates, stars_inserted, star_limit):
    """Helper function to update progress bar, abstracting HAS_TQDM checks.
    
    Args:
        pbar: Progress bar object (tqdm instance or None)
        chunk_num: Current chunk number (0-indexed)
        chunk_inserted: Number of new stars inserted in this chunk
        duplicates: Number of duplicate stars in this chunk
        stars_inserted: Total number of stars inserted so far
        star_limit: Maximum number of stars to download (None = no limit)
    """
    if HAS_TQDM and pbar is not None:
        if star_limit is not None:
            # Don't update beyond total
            remaining = max(0, star_limit - pbar.n)
            if remaining > 0:
                pbar.update(min(chunk_inserted, remaining))
        else:
            pbar.update(chunk_inserted)
        pbar.set_postfix({
            'chunk': chunk_num + 1,
            'new': f'{chunk_inserted:,}', 
            'dups': f'{duplicates:,}',
            'total': f'{stars_inserted:,}'
        })
    else:
        if star_limit is not None:
            percent = (stars_inserted / star_limit) * 100 if star_limit > 0 else 0
            dup_msg = f", {duplicates:,} duplicates" if duplicates > 0 else ""
            print(f"  ✓ Chunk {chunk_num + 1}: Inserted {chunk_inserted:,} new stars{dup_msg} (Total: {stars_inserted:,}/{star_limit:,} stars, {percent:.1f}%)")
        else:
            dup_msg = f", {duplicates:,} duplicates" if duplicates > 0 else ""
            print(f"  ✓ Chunk {chunk_num + 1}: Inserted {chunk_inserted:,} new stars{dup_msg} (Total: {stars_inserted:,} stars)")
        sys.stdout.flush()

def _download_from_gaia(star_limit=None, existing_count=0):
    """Download data from Gaia API and store in SQLite cache via shared helper."""
    from lib.gaia_client import _download_from_gaia as _gaia_download_from_gaia

    return _gaia_download_from_gaia(
        cache_db=CACHE_DB,
        chunk_size=CHUNK_SIZE,
        star_limit=star_limit,
        existing_count=existing_count,
    )


def _download_from_vizier(star_limit=None):
    """Download a roughly equivalent amount of Gaia DR3 data from VizieR and merge into cache.

    Uses the same star_limit as the Gaia request (or DEFAULT_STAR_LIMIT if None)
    so the cache is filled from both sources. Rows are merged with INSERT OR IGNORE.
    Returns number of new rows inserted from VizieR.
    """
    from lib.vizier_client import GAIA_VIZIER_CATALOG, download_vizier_catalog

    vizier_limit = star_limit if star_limit is not None else DEFAULT_STAR_LIMIT
    print(f"\nDownloading ~{vizier_limit:,} stars from VizieR ({GAIA_VIZIER_CATALOG})...")
    try:
        n = download_vizier_catalog(
            GAIA_VIZIER_CATALOG,
            CACHE_DB,
            row_limit=vizier_limit,
            merge_into_gaia=True,
        )
        print(f"  VizieR: inserted {n:,} new rows into cache.")
        return n
    except Exception as e:
        print(f"  VizieR download failed (continuing with Gaia cache only): {e}", file=sys.stderr)
        return 0

def process_star_chunk(chunk_data, target_3d, dump_positions=False, magnitude_limit=None):
    """Process a chunk of stars and return valid stars for plotting or dumping.
    
    Handles both nearby stars (with valid parallax) and background stars
    (with 0 or invalid parallax) by using a large assumed distance.
    
    If dump_positions=True, magnitude filter is skipped so all valid stars
    are returned for database output.
    
    magnitude_limit: faintest magnitude to include (uses VISIBLE_MAG_LIMIT if None).
    """
    # Extract values and handle units properly
    source_ids = np.asarray(chunk_data['source_id'])
    parallax_values = chunk_data['parallax'].value if hasattr(chunk_data['parallax'], 'value') else chunk_data['parallax']
    ra_values = chunk_data['ra'].value if hasattr(chunk_data['ra'], 'value') else chunk_data['ra']
    dec_values = chunk_data['dec'].value if hasattr(chunk_data['dec'], 'value') else chunk_data['dec']
    mag_values = chunk_data['phot_g_mean_mag'].value if hasattr(chunk_data['phot_g_mean_mag'], 'value') else chunk_data['phot_g_mean_mag']
    bp_rp_values = chunk_data['bp_rp'] if 'bp_rp' in chunk_data.colnames else np.full(len(source_ids), np.nan)
    
    # Build master mask: combine all initial filtering conditions to minimize memory copies
    # Filter for valid RA/Dec and magnitude (must be finite, but no magnitude limit)
    master_mask = (
        np.isfinite(ra_values) & np.isfinite(dec_values) & 
        np.isfinite(mag_values)  # Only check that magnitude is finite, no upper limit
    )
    
    # Handle parallax: use actual distance if valid, otherwise use large assumed distance for background stars
    # Distance (pc) = 1000 / parallax (mas) by definition of parsec; Gaia gives parallax in mas.
    # Background stars (parallax <= 0 or invalid) are very far away - use ASSUMED_BACKGROUND_DISTANCE_PC
    valid_parallax_mask = (parallax_values > 0) & np.isfinite(parallax_values)
    d_earth_pc = np.where(
        valid_parallax_mask,
        1000.0 / parallax_values,  # d (pc) = 1000 / parallax (mas)
        ASSUMED_BACKGROUND_DISTANCE_PC
    )
    
    # Add distance validation to master mask
    master_mask = master_mask & (d_earth_pc > 0) & np.isfinite(d_earth_pc) & (d_earth_pc < MAX_DISTANCE_PC)
    
    if not np.any(master_mask):
        return None
    
    # Apply master mask once to get filtered arrays for calculations
    source_ids = source_ids[master_mask]
    d_earth_pc = d_earth_pc[master_mask]
    ra_values = ra_values[master_mask]
    dec_values = dec_values[master_mask]
    mag_values = mag_values[master_mask]
    bp_rp_values = bp_rp_values[master_mask]
    has_valid_parallax = valid_parallax_mask[master_mask]
    
    # Identify background objects: stars without valid parallax OR stars beyond distance threshold
    # These should be treated as fixed background points (same position regardless of target)
    background_mask = ~has_valid_parallax | (d_earth_pc >= BACKGROUND_DISTANCE_THRESHOLD_PC)
    bright_background_mask = np.zeros(len(background_mask), dtype=bool)
    if np.any(background_mask):
        # Calculate absolute magnitude assuming the background distance
        # M = m - 5*log10(d) + 5, where d is the assumed background distance
        background_d = ASSUMED_BACKGROUND_DISTANCE_PC
        M_background = mag_values[background_mask] - 5 * np.log10(background_d) + 5
        # Mark as bright background if absolute magnitude <= threshold
        bright_background_mask[background_mask] = M_background <= BRIGHT_BACKGROUND_ABS_MAG_THRESHOLD
    
    # For bright background objects: use Earth-centric coordinates directly (fixed background)
    # For all other stars: calculate position relative to target
    if np.any(bright_background_mask):
        # Create Earth-centric coordinates for bright background objects
        # Use unit distance (1 pc) to get direction vectors, then normalize
        bright_bg_icrs = SkyCoord(ra=ra_values[bright_background_mask]*u.deg, 
                                  dec=dec_values[bright_background_mask]*u.deg, 
                                  distance=1.0*u.pc,  # Unit distance for direction
                                  frame='icrs')
        # Get unit direction vectors (normalized)
        bright_bg_cart = bright_bg_icrs.cartesian.xyz.value.T  # Shape: (N, 3)
        # Normalize to ensure unit length (should already be ~1, but ensure it)
        bright_bg_norms = np.linalg.norm(bright_bg_cart, axis=1, keepdims=True)
        bright_bg_cart = bright_bg_cart / np.where(bright_bg_norms > 1e-10, bright_bg_norms, 1.0)  # Unit direction vectors
    
    # Create 3D vectors relative to Earth in ICRS frame for all stars
    stars_earth_icrs = SkyCoord(ra=ra_values*u.deg, dec=dec_values*u.deg, 
                                distance=d_earth_pc*u.pc, frame='icrs')
    
    # Get cartesian coordinates in ICRS frame (both stars and target)
    stars_cart = stars_earth_icrs.cartesian.xyz.value.T  # Shape: (N, 3) - positions from Earth
    target_cart = target_3d.cartesian.xyz.value  # Shape: (3,) - target position from Earth
    
    # Calculate vectors from target to each star (in ICRS cartesian coordinates)
    # For bright background objects, use their Earth-centric direction directly
    # (they appear at the same projected position regardless of target star)
    vectors_from_target = stars_cart - target_cart  # Shape: (N, 3) - NumPy broadcasts (3,) to (N, 3)
    if np.any(bright_background_mask):
        # Replace with Earth-centric unit vectors for bright background objects
        # These appear at the same position regardless of target star (fixed background)
        # The azimuth/elevation calculated from these will be Earth-centric
        vectors_from_target[bright_background_mask] = bright_bg_cart
    
    # Calculate distance from target star
    # For bright background objects, distance is effectively infinite (use large value)
    d_new_pc = np.linalg.norm(vectors_from_target, axis=1)  # Euclidean distance
    if np.any(bright_background_mask):
        # Bright background objects are at effectively infinite distance
        d_new_pc[bright_background_mask] = ASSUMED_BACKGROUND_DISTANCE_PC
    
    # Build second master mask for target distance validation
    # (we need to do calculations first, then filter)
    master_mask_2 = (d_new_pc > 0) & np.isfinite(d_new_pc) & (d_new_pc < MAX_DISTANCE_PC)
    if not np.any(master_mask_2):
        return None
    
    # Apply second master mask once
    source_ids = source_ids[master_mask_2]
    d_new_pc = d_new_pc[master_mask_2]
    vectors_from_target = vectors_from_target[master_mask_2]
    ra_values = ra_values[master_mask_2]
    dec_values = dec_values[master_mask_2]
    mag_values = mag_values[master_mask_2]
    bp_rp_values = bp_rp_values[master_mask_2]
    d_earth_pc = d_earth_pc[master_mask_2]
    has_valid_parallax = has_valid_parallax[master_mask_2]
    bright_background_mask = bright_background_mask[master_mask_2]
    # Keep background_mask in sync with other arrays to avoid index mismatches
    background_mask = background_mask[master_mask_2]
    
    # Calculate angular coordinates directly from cartesian vectors
    # Use z-component to determine hemisphere: positive z = north, negative z = south
    # For plotting: azimuth from x,y components, elevation from z component
    
    # Normalize vectors to get unit direction vectors
    vector_norms = np.linalg.norm(vectors_from_target, axis=1)
    # Avoid division by zero
    vector_norms = np.where(vector_norms > 1e-10, vector_norms, 1.0)
    unit_vectors = vectors_from_target / vector_norms[:, np.newaxis]
    
    # Extract components
    x = unit_vectors[:, 0]
    y = unit_vectors[:, 1]
    z = unit_vectors[:, 2]
    
    # Calculate azimuth (angle around, 0 to 2π) from x, y components
    azimuth_rad = np.arctan2(y, x)  # Returns [-π, π]
    azimuth_rad = np.mod(azimuth_rad, 2 * np.pi)  # Normalize to [0, 2π)
    
    # Calculate elevation (angle from xy-plane, -π/2 to π/2)
    # elevation = arcsin(z) where z is the normalized z-component
    elevation_rad = np.arcsin(np.clip(z, -1.0, 1.0))  # Clip to avoid numerical errors
    
    # Calculate apparent magnitude from target star's perspective
    # Centralized magnitude logic: all stars get proper distance modulus calculation
    # Formula: m_new = M_intrinsic + 5*log10(d_new) - 5
    # where M_intrinsic = m_earth - 5*log10(d_earth) + 5
    
    # For all stars (including background), calculate absolute magnitude first
    # Then apply distance modulus for the new observer position (target star)
    
    # Calculate absolute magnitude for all stars
    # For stars with valid parallax: use actual distance from Earth
    # For background stars: use assumed background distance
    M_intrinsic = mag_values - 5 * np.log10(d_earth_pc) + 5
    
    # Calculate apparent magnitude from target star's perspective
    # This applies the distance modulus relative to the new observer position
    m_new = M_intrinsic + 5 * np.log10(d_new_pc) - 5
    
    # Special case: ALL background objects (not just bright ones) are treated as fixed background
    # They appear at the same position regardless of target, so use Earth's apparent magnitude
    # This prevents background stars from being incorrectly dimmed by distance calculations
    background_idx = np.where(background_mask)[0]
    if len(background_idx) > 0:
        m_new[background_idx] = mag_values[background_idx]
    
    # Build final master mask: combine coordinate validation, magnitude filter, and final validation
    # Verify we got valid coordinates
    valid_coord_check = (
        np.isfinite(azimuth_rad) & np.isfinite(elevation_rad) &
        np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    )
    
    # Filter by apparent magnitude (visible stars) and valid coordinates, unless dumping
    mag_limit = magnitude_limit if magnitude_limit is not None else VISIBLE_MAG_LIMIT
    if dump_positions:
        final_mask = np.isfinite(m_new) & valid_coord_check
    else:
        final_mask = (m_new < mag_limit) & np.isfinite(m_new) & valid_coord_check
    
    if not np.any(final_mask):
        return None
    
    # Apply final mask once to extract all values
    # Extract coordinates (x,y,z in pc from vectors_from_target)
    x_pc = vectors_from_target[final_mask, 0]
    y_pc = vectors_from_target[final_mask, 1]
    z_pc = vectors_from_target[final_mask, 2]
    azimuth_rad_final = azimuth_rad[final_mask]
    elevation_rad_final = elevation_rad[final_mask]
    z_final = z[final_mask]  # unit z-component for hemisphere determination
    m_plot = m_new[final_mask]
    bp_rp_final = bp_rp_values[final_mask]
    source_ids_final = source_ids[final_mask]
    
    # Final validation: ensure all extracted values are finite (combine into single check)
    valid_final = (
        np.isfinite(azimuth_rad_final) & np.isfinite(elevation_rad_final) &
        np.isfinite(z_final) & np.isfinite(m_plot) &
        np.isfinite(x_pc) & np.isfinite(y_pc) & np.isfinite(z_pc)
    )
    if not np.any(valid_final):
        return None
    
    # Distance from target (for draw order: farthest first)
    d_new_pc_final = d_new_pc[final_mask][valid_final]
    
    # Apply final validation mask once
    return {
        'source_id': source_ids_final[valid_final],
        'x_pc': x_pc[valid_final],
        'y_pc': y_pc[valid_final],
        'z_pc': z_pc[valid_final],
        'azimuth_rad': azimuth_rad_final[valid_final],
        'elevation_rad': elevation_rad_final[valid_final],
        'z': z_final[valid_final],  # unit z: positive = north, negative = south
        'm_new': m_plot[valid_final],
        'bp_rp': bp_rp_final[valid_final],
        'distance_pc': d_new_pc_final,
    }

def get_bright_galaxies():
    """Return static bright galaxy data from shared star_data module."""
    return _stardata_get_bright_galaxies()

def get_magellanic_clouds():
    """Return Magellanic Clouds data from shared star_data module."""
    return _stardata_get_magellanic_clouds()

def get_bright_stars():
    """Return bright star data from shared star_data module."""
    return _stardata_get_bright_stars()

def get_bright_deep_sky_objects():
    """Return DSO data from shared star_data module."""
    return _stardata_get_bright_deep_sky_objects()

def plot_galaxy_on_hemisphere(ax, galaxy, azimuth_rad, elevation_rad, is_north_hemisphere, clip_to_equator=False):
    """Plot a galaxy as an ellipse on a polar plot hemisphere.
    
    Args:
        ax: matplotlib polar axes
        galaxy: dict with galaxy data (from get_bright_galaxies())
        azimuth_rad: azimuth of galaxy center in radians
        elevation_rad: elevation of galaxy center in radians (may be adjusted for clipping)
        is_north_hemisphere: True for north, False for south
        clip_to_equator: If True, only plot the portion in this hemisphere
    """
    # Calculate radial position on polar plot
    if is_north_hemisphere:
        radial_center = 0.5 * np.pi - elevation_rad  # Pole at center
        equator_radial = 0.5 * np.pi  # Equator is at this radial value
    else:
        radial_center = 0.5 * np.pi + elevation_rad  # Pole at center, elevation is negative
        equator_radial = 0.5 * np.pi  # Equator is at this radial value
    
    # Convert angular size to approximate size in plot units
    # For small angles, angular size in degrees ≈ size in radial units (radians)
    major_size_rad = np.radians(galaxy['major_axis_deg'])
    minor_size_rad = np.radians(galaxy['minor_axis_deg'])
    
    # Create ellipse by plotting points around perimeter
    # Use parametric ellipse equation in local coordinates
    n_points = 256  # More points for smoother ellipse
    t = np.linspace(0, 2*np.pi, n_points)
    
    # Ellipse in local coordinates (semi-major and semi-minor axes)
    a = major_size_rad / 2  # Semi-major axis
    b = minor_size_rad / 2  # Semi-minor axis
    
    # Parametric ellipse
    x_local = a * np.cos(t)
    y_local = b * np.sin(t)
    
    # Rotate by position angle
    pa_rad = np.radians(galaxy['position_angle_deg'])
    cos_pa = np.cos(pa_rad)
    sin_pa = np.sin(pa_rad)
    x_rot = x_local * cos_pa - y_local * sin_pa
    y_rot = x_local * sin_pa + y_local * cos_pa
    
    # Convert to polar plot coordinates using proper spherical trigonometry
    # x_rot and y_rot are angular offsets in radians in the local tangent plane
    # Center point: (azimuth_rad, elevation_rad)
    # Use spherical trigonometry to compute new (azimuth, elevation) for each point
    
    # Convert center to unit vector on sphere
    center_x = np.cos(elevation_rad) * np.cos(azimuth_rad)
    center_y = np.cos(elevation_rad) * np.sin(azimuth_rad)
    center_z = np.sin(elevation_rad)
    center_vec = np.array([center_x, center_y, center_z])
    
    # Create local tangent frame vectors at center
    # North vector (increasing elevation): derivative w.r.t. elevation
    north_x = -np.sin(elevation_rad) * np.cos(azimuth_rad)
    north_y = -np.sin(elevation_rad) * np.sin(azimuth_rad)
    north_z = np.cos(elevation_rad)
    north_vec = np.array([north_x, north_y, north_z])
    
    # East vector (increasing azimuth): derivative w.r.t. azimuth, scaled by cos(elevation)
    east_x = -np.sin(azimuth_rad)
    east_y = np.cos(azimuth_rad)
    east_z = 0.0
    east_vec = np.array([east_x, east_y, east_z])
    
    # For each ellipse point, compute offset vector in tangent plane
    # x_rot is offset in north direction, y_rot is offset in east direction
    # Use small-angle approximation for tangent plane (valid for small angular sizes)
    # Then project back to unit sphere
    
    # Compute offset vectors for all points at once
    offset_vecs = (x_rot[:, np.newaxis] * north_vec + 
                   y_rot[:, np.newaxis] * east_vec)
    
    # Project to unit sphere: new_vec = normalize(center_vec + offset_vec)
    new_vecs = center_vec + offset_vecs
    new_norms = np.linalg.norm(new_vecs, axis=1, keepdims=True)
    new_norms = np.where(new_norms > 1e-10, new_norms, 1.0)  # Avoid division by zero
    new_vecs_normalized = new_vecs / new_norms
    
    # Convert back to spherical coordinates (azimuth, elevation)
    new_x = new_vecs_normalized[:, 0]
    new_y = new_vecs_normalized[:, 1]
    new_z = new_vecs_normalized[:, 2]
    
    new_azimuth = np.arctan2(new_y, new_x)
    new_azimuth = np.mod(new_azimuth, 2 * np.pi)  # Wrap to [0, 2π)
    new_elevation = np.arcsin(np.clip(new_z, -1.0, 1.0))
    
    # Convert to polar plot coordinates
    if is_north_hemisphere:
        ellipse_radial = 0.5 * np.pi - new_elevation
    else:
        ellipse_radial = 0.5 * np.pi + new_elevation
    
    ellipse_azimuth = new_azimuth
    
    # Filter points within valid hemisphere bounds
    if is_north_hemisphere:
        if clip_to_equator:
            # Only show points above equator (radial <= equator_radial for north)
            valid = ellipse_radial <= equator_radial
        else:
            valid = ellipse_radial >= 0
    else:
        if clip_to_equator:
            # Only show points below equator (radial >= equator_radial for south)
            valid = ellipse_radial >= equator_radial
        else:
            valid = ellipse_radial >= 0
    
    if np.any(valid):
        # If clipping, add equator intersection points
        if clip_to_equator and np.sum(valid) < len(valid):
            # Find where ellipse crosses equator and add intersection points
            # This is a simplified approach - add points at equator boundary
            valid_indices = np.where(valid)[0]
            if len(valid_indices) > 0:
                # Add points at boundaries
                boundary_points_az = []
                boundary_points_rad = []
                
                # Check transitions from invalid to valid and vice versa
                for i in range(len(valid)):
                    if valid[i] and not valid[(i-1) % len(valid)]:
                        # Entering valid region - add point at equator
                        boundary_points_az.append(ellipse_azimuth[i])
                        boundary_points_rad.append(equator_radial)
                    elif not valid[i] and valid[(i-1) % len(valid)]:
                        # Leaving valid region - add point at equator
                        boundary_points_az.append(ellipse_azimuth[(i-1) % len(valid)])
                        boundary_points_rad.append(equator_radial)
                
                # Combine valid points with boundary points
                if boundary_points_az:
                    all_az = np.append(ellipse_azimuth[valid], boundary_points_az)
                    all_rad = np.append(ellipse_radial[valid], boundary_points_rad)
                    # Sort by angle to maintain order
                    sort_idx = np.argsort(all_az)
                    ellipse_azimuth_plot = all_az[sort_idx]
                    ellipse_radial_plot = all_rad[sort_idx]
                else:
                    ellipse_azimuth_plot = ellipse_azimuth[valid]
                    ellipse_radial_plot = ellipse_radial[valid]
            else:
                ellipse_azimuth_plot = ellipse_azimuth[valid]
                ellipse_radial_plot = ellipse_radial[valid]
        else:
            ellipse_azimuth_plot = ellipse_azimuth[valid]
            ellipse_radial_plot = ellipse_radial[valid]
        
        # Plot ellipse outline
        ax.plot(ellipse_azimuth_plot, ellipse_radial_plot, 
                color='cyan', linewidth=2, alpha=0.8, transform=ax.transData)
        
        # Add label at center (only if center is in this hemisphere)
        # Offset down for cyan galaxies
        center_in_hemisphere = (is_north_hemisphere and elevation_rad > 0) or (not is_north_hemisphere and elevation_rad < 0)
        if center_in_hemisphere or not clip_to_equator:
            fontsize = _label_font_size(LABEL_FONT_SIZE_RATIO_LARGE)
            ann = ax.annotate(galaxy['name'],
                             xy=(azimuth_rad, radial_center),
                             xytext=(0, -10),  # 10 points below the star
                             textcoords='offset points',
                             ha='center', va='top',
                             color='cyan', fontsize=fontsize, weight='bold',
                             transform=ax.transData)
            ann.set_path_effects(_text_stroke_effects())

def get_relative_coords(target_xyz, object_xyz):
    """Delegate to shared math helper."""
    return _math_get_relative_coords(target_xyz, object_xyz)


def transform_to_target_frame(target_coord, object_ra_deg, object_dec_deg, object_distance_pc, use_earth_centric_approx=False):
    """Delegate to shared math helper."""
    return _math_transform_to_target_frame(
        target_coord,
        object_ra_deg,
        object_dec_deg,
        object_distance_pc,
        use_earth_centric_approx=use_earth_centric_approx,
    )

def calculate_galaxy_coordinates(galaxy, target_3d):
    """Calculate azimuth and elevation of a galaxy from target star's perspective.
    
    Uses distance estimates from galaxy data.
    
    Returns: (azimuth_rad, elevation_rad, z)
    """
    # Use distance from galaxy data if available, otherwise use large assumed distance
    galaxy_distance_pc = galaxy.get('distance_pc', 100000.0)
    
    return transform_to_target_frame(
        target_3d,
        galaxy['ra_deg'],
        galaxy['dec_deg'],
        galaxy_distance_pc,
        use_earth_centric_approx=False
    )

def calculate_magellanic_cloud_coordinates(mc, target_3d):
    """Calculate azimuth, elevation, and apparent angular size of a Magellanic Cloud from target star's perspective.
    
    The apparent size changes based on the distance from the target star to the galaxy.
    Uses physical dimensions to calculate apparent angular size.
    
    Args:
        mc: Magellanic Cloud dict with physical_major_axis_pc, physical_minor_axis_pc, distance_pc
        target_3d: SkyCoord of target star position
    
    Returns: (azimuth_rad, elevation_rad, z, apparent_major_axis_deg, apparent_minor_axis_deg)
    """
    # Get Earth-centric coordinates
    mc_earth_icrs = SkyCoord(ra=mc['ra_deg']*u.deg, dec=mc['dec_deg']*u.deg, 
                             distance=mc['distance_pc']*u.pc, frame='icrs')
    mc_earth_cart = mc_earth_icrs.cartesian.xyz.value  # Shape: (3,)
    
    # Get target position
    target_cart = target_3d.cartesian.xyz.value  # Shape: (3,)
    
    # Calculate distance from target to galaxy
    vector_from_target = mc_earth_cart - target_cart
    distance_from_target_pc = np.linalg.norm(vector_from_target)
    
    # Safety check: avoid division by zero or extremely small distances
    # If target is very close to galaxy, use Earth-centric angular size as fallback
    min_distance_pc = 1.0  # Minimum distance threshold (1 pc)
    if distance_from_target_pc < min_distance_pc:
        distance_from_target_pc = mc['distance_pc']  # Use Earth distance as fallback
        # Recalculate vector using Earth-centric distance
        mc_earth_icrs_unit = SkyCoord(ra=mc['ra_deg']*u.deg, dec=mc['dec_deg']*u.deg, 
                                      distance=1.0*u.pc, frame='icrs')
        mc_earth_unit_cart = mc_earth_icrs_unit.cartesian.xyz.value
        vector_from_target = mc_earth_unit_cart * mc['distance_pc'] - target_cart
        distance_from_target_pc = np.linalg.norm(vector_from_target)
        if distance_from_target_pc < min_distance_pc:
            distance_from_target_pc = mc['distance_pc']
    
    # Calculate apparent angular size from target star's perspective
    # Angular size (radians) = physical size (pc) / distance (pc)
    # Convert to degrees
    apparent_major_axis_rad = mc['physical_major_axis_pc'] / distance_from_target_pc
    apparent_minor_axis_rad = mc['physical_minor_axis_pc'] / distance_from_target_pc
    apparent_major_axis_deg = np.degrees(apparent_major_axis_rad)
    apparent_minor_axis_deg = np.degrees(apparent_minor_axis_rad)
    
    # Calculate direction (azimuth, elevation) from target to galaxy
    unit_vector = vector_from_target / distance_from_target_pc
    x, y, z = unit_vector
    
    azimuth_rad = np.arctan2(y, x)
    azimuth_rad = np.mod(azimuth_rad, 2 * np.pi)  # Normalize to [0, 2π)
    elevation_rad = np.arcsin(np.clip(z, -1.0, 1.0))
    
    return azimuth_rad, elevation_rad, z, apparent_major_axis_deg, apparent_minor_axis_deg

def plot_magellanic_cloud(ax, mc, azimuth_rad, elevation_rad, apparent_major_axis_deg, 
                          apparent_minor_axis_deg, is_north_hemisphere, clip_to_equator=False):
    """Plot a Magellanic Cloud as an ellipse on a polar plot hemisphere with dynamic sizing.
    
    Args:
        ax: matplotlib polar axes
        mc: dict with Magellanic Cloud data
        azimuth_rad: azimuth of cloud center in radians
        elevation_rad: elevation of cloud center in radians
        apparent_major_axis_deg: apparent major axis in degrees (from target star's perspective)
        apparent_minor_axis_deg: apparent minor axis in degrees (from target star's perspective)
        is_north_hemisphere: True for north, False for south
        clip_to_equator: If True, only plot the portion in this hemisphere
    """
    # Use the same ellipse plotting logic as plot_galaxy_on_hemisphere, but with dynamic sizes
    # Calculate radial position on polar plot
    if is_north_hemisphere:
        radial_center = 0.5 * np.pi - elevation_rad  # Pole at center
        equator_radial = 0.5 * np.pi  # Equator is at this radial value
    else:
        radial_center = 0.5 * np.pi + elevation_rad  # Pole at center, elevation is negative
        equator_radial = 0.5 * np.pi  # Equator is at this radial value
    
    # Convert angular size to approximate size in plot units
    major_size_rad = np.radians(apparent_major_axis_deg)
    minor_size_rad = np.radians(apparent_minor_axis_deg)
    
    # Create ellipse by plotting points around perimeter
    n_points = 256  # More points for smoother ellipse
    t = np.linspace(0, 2*np.pi, n_points)
    
    # Ellipse in local coordinates (semi-major and semi-minor axes)
    a = major_size_rad / 2  # Semi-major axis
    b = minor_size_rad / 2  # Semi-minor axis
    
    # Parametric ellipse
    x_local = a * np.cos(t)
    y_local = b * np.sin(t)
    
    # Rotate by position angle
    pa_rad = np.radians(mc['position_angle_deg'])
    cos_pa = np.cos(pa_rad)
    sin_pa = np.sin(pa_rad)
    x_rot = x_local * cos_pa - y_local * sin_pa
    y_rot = x_local * sin_pa + y_local * cos_pa
    
    # Convert to polar plot coordinates using proper spherical trigonometry
    # Convert center to unit vector on sphere
    center_x = np.cos(elevation_rad) * np.cos(azimuth_rad)
    center_y = np.cos(elevation_rad) * np.sin(azimuth_rad)
    center_z = np.sin(elevation_rad)
    center_vec = np.array([center_x, center_y, center_z])
    
    # Create local tangent frame vectors at center
    north_x = -np.sin(elevation_rad) * np.cos(azimuth_rad)
    north_y = -np.sin(elevation_rad) * np.sin(azimuth_rad)
    north_z = np.cos(elevation_rad)
    north_vec = np.array([north_x, north_y, north_z])
    
    east_x = -np.sin(azimuth_rad)
    east_y = np.cos(azimuth_rad)
    east_z = 0.0
    east_vec = np.array([east_x, east_y, east_z])
    
    # Compute offset vectors for all points at once
    offset_vecs = (x_rot[:, np.newaxis] * north_vec + 
                   y_rot[:, np.newaxis] * east_vec)
    
    # Project to unit sphere: new_vec = normalize(center_vec + offset_vec)
    new_vecs = center_vec + offset_vecs
    new_norms = np.linalg.norm(new_vecs, axis=1, keepdims=True)
    new_norms = np.where(new_norms > 1e-10, new_norms, 1.0)  # Avoid division by zero
    new_vecs_normalized = new_vecs / new_norms
    
    # Convert back to spherical coordinates (azimuth, elevation)
    new_x = new_vecs_normalized[:, 0]
    new_y = new_vecs_normalized[:, 1]
    new_z = new_vecs_normalized[:, 2]
    
    new_azimuth = np.arctan2(new_y, new_x)
    new_azimuth = np.mod(new_azimuth, 2 * np.pi)  # Wrap to [0, 2π)
    new_elevation = np.arcsin(np.clip(new_z, -1.0, 1.0))
    
    # Convert to polar plot coordinates
    if is_north_hemisphere:
        ellipse_radial = 0.5 * np.pi - new_elevation
    else:
        ellipse_radial = 0.5 * np.pi + new_elevation
    
    ellipse_azimuth = new_azimuth
    
    # Filter points within valid hemisphere bounds
    if is_north_hemisphere:
        if clip_to_equator:
            valid = ellipse_radial <= equator_radial
        else:
            valid = ellipse_radial >= 0
    else:
        if clip_to_equator:
            valid = ellipse_radial >= equator_radial
        else:
            valid = ellipse_radial >= 0
    
    if np.any(valid):
        ellipse_azimuth_plot = ellipse_azimuth[valid]
        ellipse_radial_plot = ellipse_radial[valid]
        
        # Plot ellipse outline (use cyan like other galaxies, or distinct color)
        ax.plot(ellipse_azimuth_plot, ellipse_radial_plot, 
                color='cyan', linewidth=2, alpha=0.8, transform=ax.transData)
        
        # Add label at center (only if center is in this hemisphere)
        center_in_hemisphere = (is_north_hemisphere and elevation_rad > 0) or (not is_north_hemisphere and elevation_rad < 0)
        if center_in_hemisphere or not clip_to_equator:
            fontsize = _label_font_size(LABEL_FONT_SIZE_RATIO_LARGE)
            ann = ax.annotate(mc['name'],
                             xy=(azimuth_rad, radial_center),
                             xytext=(0, -10),  # 10 points below
                             textcoords='offset points',
                             ha='center', va='top',
                             color='cyan', fontsize=fontsize, weight='bold',
                             transform=ax.transData)
            ann.set_path_effects(_text_stroke_effects())

def calculate_star_coordinates(star, target_3d):
    """Calculate azimuth and elevation of a bright star from target star's perspective.
    
    Uses parallax values when available, otherwise assumes large distance.
    For very distant stars, the direction is nearly the same as from Earth.
    
    Returns: (azimuth_rad, elevation_rad, z)
    """
    # Use parallax if available, otherwise use large assumed distance
    parallax_mas = star.get('parallax_mas', None)
    if parallax_mas and parallax_mas > 0:
        # Distance from parallax: d (pc) = 1000 / parallax (mas)
        star_distance_pc = 1000.0 / parallax_mas
    else:
        # Very distant stars without parallax measurement
        star_distance_pc = 100000.0
    
    return transform_to_target_frame(
        target_3d,
        star['ra_deg'],
        star['dec_deg'],
        star_distance_pc,
        use_earth_centric_approx=True  # Enable special handling for very distant stars
    )

def calculate_zone_of_avoidance_polygons(target_3d, disk_thickness_pc=2000.0, disk_radius_pc=50000.0, n_azimuth=180, n_elevation=90):
    """Calculate polygonal regions for the zone of avoidance (galactic dust/gas obscuration).
    
    Models the zone of avoidance as a 3D disk of gas and dust (the Milky Way disk).
    For each direction from the target star, checks if the line of sight intersects
    the galactic disk volume.
    
    Args:
        target_3d: SkyCoord of target star position (in ICRS frame)
        disk_thickness_pc: Thickness of galactic disk in parsecs (default 2000 pc)
        disk_radius_pc: Radius of galactic disk in parsecs (default 50000 pc)
        n_azimuth: Number of azimuth samples (default 360)
        n_elevation: Number of elevation samples (default 180)
    
    Returns:
        Dict with 'azimuth_rad', 'elevation_rad', 'z' arrays for directions that
        intersect the galactic disk
    """
    # Get target star's position in galactic coordinates
    target_galactic = target_3d.transform_to('galactic')
    target_cart_gal = target_galactic.cartesian.xyz.value  # Shape: (3,)
    
    # In galactic coordinates:
    # x_gal: toward galactic center
    # y_gal: direction of galactic rotation
    # z_gal: perpendicular to galactic plane (north = positive)
    
    # The disk is centered on z_gal = 0, with thickness disk_thickness_pc
    # and extends to radius disk_radius_pc in the x-y plane
    
    # Create a grid of directions from the target star
    azimuth_samples = np.linspace(0, 2*np.pi, n_azimuth, endpoint=False)
    elevation_samples = np.linspace(-np.pi/2, np.pi/2, n_elevation)
    
    # Create meshgrid
    AZ, EL = np.meshgrid(azimuth_samples, elevation_samples)
    
    # Convert to unit direction vectors in ICRS frame
    x_dir = np.cos(EL) * np.cos(AZ)
    y_dir = np.cos(EL) * np.sin(AZ)
    z_dir = np.sin(EL)
    
    # Flatten for processing
    x_dir_flat = x_dir.flatten()
    y_dir_flat = y_dir.flatten()
    z_dir_flat = z_dir.flatten()
    
    # Compute rotation matrix from ICRS to Galactic frame
    # Transform ICRS basis vectors (x, y, z) to Galactic frame to get rotation matrix
    # Use unit distance to get pure direction vectors
    unit_distance = 1.0 * u.pc
    icrs_basis = np.eye(3)  # ICRS basis vectors: [1,0,0], [0,1,0], [0,0,1]
    
    # Transform each basis vector to Galactic frame
    galactic_basis = np.zeros((3, 3))
    for i in range(3):
        # Create SkyCoord for basis vector in ICRS
        basis_icrs = SkyCoord(CartesianRepresentation(
            icrs_basis[i, 0] * unit_distance,
            icrs_basis[i, 1] * unit_distance,
            icrs_basis[i, 2] * unit_distance
        ), frame='icrs')
        # Transform to Galactic
        basis_gal = basis_icrs.transform_to('galactic')
        galactic_basis[i, :] = basis_gal.cartesian.xyz.value
    
    # Rotation matrix: columns are Galactic frame basis vectors expressed in ICRS
    # To transform ICRS -> Galactic, we use the transpose (rows are Galactic basis in ICRS)
    rotation_matrix = galactic_basis.T  # Shape: (3, 3)
    
    # Convert all direction vectors from ICRS to Galactic frame at once (vectorized)
    dir_vectors_icrs = np.column_stack([x_dir_flat, y_dir_flat, z_dir_flat])  # Shape: (N, 3)
    dir_vectors_gal = dir_vectors_icrs @ rotation_matrix  # Shape: (N, 3)
    
    # Normalize to ensure unit vectors (should already be unit, but ensure numerical precision)
    dir_norms = np.linalg.norm(dir_vectors_gal, axis=1, keepdims=True)
    dir_norms = np.where(dir_norms > 1e-10, dir_norms, 1.0)  # Avoid division by zero
    dir_vectors_gal_unit = dir_vectors_gal / dir_norms  # Shape: (N, 3)
    
    # Check intersections for all directions
    intersects_disk = []
    azimuth_result = []
    elevation_result = []
    z_result = []
    
    for i in range(len(x_dir_flat)):
        dir_gal_unit = dir_vectors_gal_unit[i, :]
        
        # Check if line of sight from target intersects the disk
        # The disk is: |z_gal| < disk_thickness_pc/2 and sqrt(x_gal^2 + y_gal^2) < disk_radius_pc
        
        # Parametric line: point = target_cart_gal + t * dir_gal_unit
        # We need to find t where the line intersects the disk
        
        # Check if direction has a z-component (perpendicular to plane)
        if abs(dir_gal_unit[2]) < 1e-10:
            # Direction is parallel to plane - check if target is in disk
            z_target = target_cart_gal[2]
            if abs(z_target) < disk_thickness_pc / 2:
                r_target = np.sqrt(target_cart_gal[0]**2 + target_cart_gal[1]**2)
                if r_target < disk_radius_pc:
                    # Target is in disk, all directions in plane intersect
                    intersects = True
                else:
                    intersects = False
            else:
                intersects = False
        else:
            # Find where line crosses z = ±disk_thickness_pc/2
            z_target = target_cart_gal[2]
            
            # Distance to top of disk (z = +disk_thickness_pc/2)
            t_top = (disk_thickness_pc/2 - z_target) / dir_gal_unit[2]
            # Distance to bottom of disk (z = -disk_thickness_pc/2)
            t_bottom = (-disk_thickness_pc/2 - z_target) / dir_gal_unit[2]
            
            # Check both intersections
            t_values = []
            if t_top > 0:
                point_top = target_cart_gal + dir_gal_unit * t_top
                r_top = np.sqrt(point_top[0]**2 + point_top[1]**2)
                if r_top < disk_radius_pc:
                    t_values.append(t_top)
            
            if t_bottom > 0:
                point_bottom = target_cart_gal + dir_gal_unit * t_bottom
                r_bottom = np.sqrt(point_bottom[0]**2 + point_bottom[1]**2)
                if r_bottom < disk_radius_pc:
                    t_values.append(t_bottom)
            
            # Also check if line intersects disk boundary (circular edge)
            # This is more complex, but for now we'll use the z-boundary check
            
            intersects = len(t_values) > 0
        
        if intersects:
            azimuth_result.append(azimuth_samples[i % n_azimuth])
            elevation_result.append(elevation_samples[i // n_azimuth])
            z_result.append(z_dir_flat[i])
    
    if len(azimuth_result) == 0:
        return {
            'azimuth_rad': np.array([]),
            'elevation_rad': np.array([]),
            'z': np.array([])
        }
    
    return {
        'azimuth_rad': np.array(azimuth_result),
        'elevation_rad': np.array(elevation_result),
        'z': np.array(z_result)
    }

def plot_zone_of_avoidance(ax, target_3d, is_north_hemisphere, disk_thickness_pc=2000.0, disk_radius_pc=50000.0):
    """Plot the zone of avoidance as dark grey shaded regions on a polar plot.
    
    Args:
        ax: matplotlib polar axes
        target_3d: SkyCoord of target star position
        is_north_hemisphere: True for north, False for south
        disk_thickness_pc: Thickness of galactic disk in parsecs
        disk_radius_pc: Radius of galactic disk in parsecs
    """
    try:
        poly = calculate_zone_of_avoidance_polygons(
            target_3d, 
            disk_thickness_pc=disk_thickness_pc,
            disk_radius_pc=disk_radius_pc
        )
        
        az = poly['azimuth_rad']
        el = poly['elevation_rad']
        z_poly = poly['z']
        
        if len(az) == 0:
            return
        
        # Filter to only include points in this hemisphere
        if is_north_hemisphere:
            mask = z_poly > 0
        else:
            mask = z_poly < 0
        
        if not np.any(mask):
            return
        
        az_filtered = az[mask]
        el_filtered = el[mask]
        
        # Convert to polar plot coordinates
        if is_north_hemisphere:
            radial = 0.5 * np.pi - el_filtered
        else:
            radial = 0.5 * np.pi + el_filtered
        
        # Create a filled region using scatter with large points or using fill
        # For better visualization, we'll use a scatter plot with large semi-transparent points
        # Or create a convex hull and fill it
        
        # Simple approach: use scatter with large points for density visualization
        # Or create a boundary polygon
        
        # Try to create a boundary by finding the convex hull or by sorting
        # For now, use a scatter-based approach with large alpha
        ax.scatter(az_filtered, radial, 
                  s=50, color='#111111', alpha=0.3, 
                  transform=ax.transData, zorder=0, edgecolors='none')
        
        # Also try to create a filled polygon by finding the boundary
        # Sort by azimuth and create a boundary
        if len(az_filtered) > 3:
            # Create a boundary polygon
            # Sort by azimuth
            sort_idx = np.argsort(az_filtered)
            az_sorted = az_filtered[sort_idx]
            radial_sorted = radial[sort_idx]
            
            # Create a simple boundary by taking min/max radial at each azimuth
            # This is a simplification, but should work for most cases
            unique_az, az_indices = np.unique(az_sorted, return_index=True)
            if len(unique_az) > 2:
                # For each unique azimuth, find min and max radial
                boundary_az = []
                boundary_rad = []
                
                for az_val in unique_az:
                    az_mask = az_sorted == az_val
                    if np.any(az_mask):
                        rad_values = radial_sorted[az_mask]
                        boundary_az.append(az_val)
                        boundary_rad.append(np.min(rad_values))
                        boundary_az.append(az_val)
                        boundary_rad.append(np.max(rad_values))
                
                if len(boundary_az) > 2:
                    # Close the polygon
                    boundary_az = np.array(boundary_az + [boundary_az[0]])
                    boundary_rad = np.array(boundary_rad + [boundary_rad[0]])
                    
                    ax.fill(boundary_az, boundary_rad,
                           color='#111111', alpha=0.3, edgecolor='#111111', linewidth=0.5,
                           transform=ax.transData, zorder=0)
    except Exception as e:
        # If zone of avoidance calculation fails, just skip it
        print(f"Warning: Could not plot zone of avoidance: {e}")
        import traceback
        traceback.print_exc()

def calculate_dso_coordinates(dso, target_3d):
    """Calculate azimuth and elevation of a deep sky object from target star's perspective.
    
    Uses distance estimates from DSO data when available.
    
    Returns: (azimuth_rad, elevation_rad, z)
    """
    # Use distance from DSO data if available, otherwise use large assumed distance
    # Most DSOs (nebulas, clusters) are at distances from hundreds to thousands of parsecs
    dso_distance_pc = dso.get('distance_pc', None)
    if dso_distance_pc is None or dso_distance_pc <= 0:
        # Default distances based on object type
        if dso.get('object_type') == 'globular':
            dso_distance_pc = 10000.0  # Globular clusters typically 5-50 kpc
        elif dso.get('object_type') == 'nebula':
            dso_distance_pc = 2000.0   # Nebulas typically 1-5 kpc
        else:  # cluster (open cluster)
            dso_distance_pc = 500.0   # Open clusters typically 100-2000 pc
    
    return transform_to_target_frame(
        target_3d,
        dso['ra_deg'],
        dso['dec_deg'],
        dso_distance_pc,
        use_earth_centric_approx=False
    )

def calculate_sol_coordinates(target_3d):
    """Calculate azimuth and elevation of Sol (Sun) from target star's perspective.
    
    Sol is at the origin (0,0,0) in ICRS coordinates, so the vector from target to Sol
    is simply -target_3d.
    
    Returns: (azimuth_rad, elevation_rad, z)
    """
    # Sol is at origin (0, 0, 0)
    sol_cart = np.array([0.0, 0.0, 0.0])
    target_cart = target_3d.cartesian.xyz.value  # Shape: (3,)
    
    # Use unified get_relative_coords function
    _, _, unit_vector = get_relative_coords(target_cart, sol_cart)
    
    # Check if target is at origin (Sol itself)
    if np.linalg.norm(unit_vector) < 1e-10:
        return 0.0, 0.0, 0.0
    
    # Convert to azimuth and elevation
    x, y, z = unit_vector
    azimuth_rad = np.arctan2(y, x)
    azimuth_rad = np.mod(azimuth_rad, 2 * np.pi)  # Normalize to [0, 2π)
    elevation_rad = np.arcsin(np.clip(z, -1.0, 1.0))
    
    return azimuth_rad, elevation_rad, z

def calculate_sagittarius_a_coordinates(target_3d):
    """Calculate azimuth and elevation of Sagittarius A* (galactic center) from target star's perspective.

    The direction to Sgr A* depends on the target: we form the vector from the target's
    3D position to Sgr A*'s 3D position (both in ICRS), then convert to azimuth/elevation
    in the same convention used for the star field (ICRS north = positive z).
    Sagittarius A* fixed position (from Earth):
    - RA: 266.4168°, Dec: -29.0078°, Distance: ~8,000 pc (~26,000 ly)

    Returns: (azimuth_rad, elevation_rad, z) as plain floats; position is target-relative.
    """
    # Sgr A* position in ICRS (Cartesian, origin at Solar System barycenter)
    sgr_a_ra_deg = 266.4168
    sgr_a_dec_deg = -29.0078
    sgr_a_distance_pc = 8000.0  # ~8 kpc from Earth
    sgr_a_icrs = SkyCoord(
        ra=sgr_a_ra_deg * u.deg,
        dec=sgr_a_dec_deg * u.deg,
        distance=sgr_a_distance_pc * u.pc,
        frame="icrs",
    )
    sgr_a_cart = np.asarray(sgr_a_icrs.cartesian.xyz.value)   # (3,) in pc
    target_cart = np.asarray(target_3d.cartesian.xyz.value)  # (3,) in pc

    # Vector from target to Sgr A* (direction in which to look from target)
    vector_from_target = sgr_a_cart - target_cart
    distance = np.linalg.norm(vector_from_target)
    if distance < 1e-10:
        return 0.0, 0.0, 0.0  # Target is at Sgr A*
    unit_vector = vector_from_target / distance
    x, y, z = float(unit_vector[0]), float(unit_vector[1]), float(unit_vector[2])

    azimuth_rad = float(np.arctan2(y, x))
    azimuth_rad = float(np.mod(azimuth_rad, 2 * np.pi))
    elevation_rad = float(np.arcsin(np.clip(z, -1.0, 1.0)))
    return azimuth_rad, elevation_rad, z

def plot_sagittarius_a_reference(ax, azimuth_rad, elevation_rad, is_north_hemisphere):
    """Plot Sagittarius A* as a special reference point on a polar plot hemisphere.
    
    Args:
        ax: matplotlib polar axes
        azimuth_rad: azimuth of Sagittarius A* in radians
        elevation_rad: elevation of Sagittarius A* in radians
        is_north_hemisphere: True for north, False for south
    """
    # Calculate radial position on polar plot
    if is_north_hemisphere:
        radial = 0.5 * np.pi - elevation_rad  # Pole at center
    else:
        radial = 0.5 * np.pi + elevation_rad  # Pole at center, elevation is negative
    
    # Only plot if in this hemisphere
    if (is_north_hemisphere and elevation_rad > 0) or (not is_north_hemisphere and elevation_rad < 0):
        # Plot Sagittarius A* as a special marker (larger, different color)
        ax.scatter(azimuth_rad, radial, s=300, color='red', 
                  alpha=1.0, edgecolors='darkred', linewidths=3, 
                  marker='*', transform=ax.transData, zorder=10)
        
        # Add label using ax.annotate with offset points
        fontsize = _label_font_size(LABEL_FONT_SIZE_RATIO_LARGE)
        ann = ax.annotate('Sgr A*',
                         xy=(azimuth_rad, radial),
                         xytext=(0, -10),  # 10 points below the star
                         textcoords='offset points',
                         ha='center', va='top',
                         color='red', fontsize=fontsize, weight='bold',
                         transform=ax.transData, zorder=10)
        ann.set_path_effects(_text_stroke_effects())

def plot_sol_reference(ax, azimuth_rad, elevation_rad, is_north_hemisphere):
    """Plot Sol (Sun) as a special reference star on a polar plot hemisphere.
    
    Args:
        ax: matplotlib polar axes
        azimuth_rad: azimuth of Sol in radians
        elevation_rad: elevation of Sol in radians
        is_north_hemisphere: True for north, False for south
    """
    # Calculate radial position on polar plot
    if is_north_hemisphere:
        radial = 0.5 * np.pi - elevation_rad  # Pole at center
    else:
        radial = 0.5 * np.pi + elevation_rad  # Pole at center, elevation is negative
    
    # Only plot if in this hemisphere
    if (is_north_hemisphere and elevation_rad > 0) or (not is_north_hemisphere and elevation_rad < 0):
        # Plot Sol as a special marker (larger, different color)
        ax.scatter(azimuth_rad, radial, s=300, color='orange', 
                  alpha=1.0, edgecolors='red', linewidths=3, 
                  marker='*', transform=ax.transData, zorder=10)
        
        # Add label using ax.annotate with offset points
        fontsize = _label_font_size(LABEL_FONT_SIZE_RATIO_LARGE)
        ann = ax.annotate('Sol',
                         xy=(azimuth_rad, radial),
                         xytext=(0, -10),  # 10 points below the star
                         textcoords='offset points',
                         ha='center', va='top',
                         color='orange', fontsize=fontsize, weight='bold',
                         transform=ax.transData, zorder=10)
        ann.set_path_effects(_text_stroke_effects())

def calculate_point_size_by_magnitude(magnitude, min_size=1, max_size=10, magnitude_brightest=10):
    """Delegate to shared point-size helper."""
    return _plot_calculate_point_size_by_magnitude(
        magnitude,
        min_size=min_size,
        max_size=max_size,
        magnitude_brightest=magnitude_brightest,
    )

def plot_star_label(ax, star, azimuth_rad, elevation_rad, is_north_hemisphere, apparent_mag_from_target, magnitude_limit=None, point_size_min=None, point_size_max=None):
    """Plot a bright star with label on a polar plot hemisphere.
    
    Args:
        ax: matplotlib polar axes
        star: dict with star data (from get_bright_stars())
        azimuth_rad: azimuth of star in radians
        elevation_rad: elevation of star in radians
        is_north_hemisphere: True for north, False for south
        apparent_mag_from_target: Apparent magnitude from target star's perspective
        magnitude_limit: faintest magnitude for sizing (uses VISIBLE_MAG_LIMIT if None)
        point_size_min: smallest point size (uses STAR_POINT_MIN_SIZE if None)
        point_size_max: largest point size (uses STAR_POINT_MAX_SIZE if None)
    """
    mag_limit = magnitude_limit if magnitude_limit is not None else VISIBLE_MAG_LIMIT
    min_pt = point_size_min if point_size_min is not None else STAR_POINT_MIN_SIZE
    max_pt = point_size_max if point_size_max is not None else STAR_POINT_MAX_SIZE
    # Calculate radial position on polar plot
    if is_north_hemisphere:
        radial = 0.5 * np.pi - elevation_rad  # Pole at center
    else:
        radial = 0.5 * np.pi + elevation_rad  # Pole at center, elevation is negative
    
    # Only plot if in this hemisphere
    if (is_north_hemisphere and elevation_rad > 0) or (not is_north_hemisphere and elevation_rad < 0):
        # Calculate point size based on apparent magnitude from target
        # Use same size range as background stars for consistency
        point_size = calculate_point_size_by_magnitude(
            apparent_mag_from_target,
            min_size=min_pt,
            max_size=max_pt,
            magnitude_brightest=mag_limit
        )
        
        # Only plot if point size is greater than 0
        if point_size > 0:
            # Plot star as a point (rendered on top of background stars)
            ax.scatter(azimuth_rad, radial, s=point_size, color='yellow', 
                      alpha=0.9, edgecolors='orange', linewidths=1, transform=ax.transData, zorder=10)
        
        # Add label using ax.annotate with offset points
        fontsize = _label_font_size(LABEL_FONT_SIZE_RATIO_SMALL)
        ann = ax.annotate(star['name'],
                         xy=(azimuth_rad, radial),
                         xytext=(0, -10),  # 10 points below the star
                         textcoords='offset points',
                         ha='center', va='top',
                         color='yellow', fontsize=fontsize, weight='bold',
                         transform=ax.transData)
        ann.set_path_effects(_text_stroke_effects())

def plot_dso_on_hemisphere(ax, dso, azimuth_rad, elevation_rad, is_north_hemisphere, clip_to_equator=False):
    """Plot a deep sky object (nebula, cluster) as an ellipse on a polar plot hemisphere.
    
    Args:
        ax: matplotlib polar axes
        dso: dict with DSO data (from get_bright_deep_sky_objects())
        azimuth_rad: azimuth of DSO center in radians
        elevation_rad: elevation of DSO center in radians
        is_north_hemisphere: True for north, False for south
        clip_to_equator: If True, only plot the portion in this hemisphere
    """
    # Calculate radial position on polar plot
    if is_north_hemisphere:
        radial_center = 0.5 * np.pi - elevation_rad  # Pole at center
        equator_radial = 0.5 * np.pi  # Equator is at this radial value
    else:
        radial_center = 0.5 * np.pi + elevation_rad  # Pole at center, elevation is negative
        equator_radial = 0.5 * np.pi  # Equator is at this radial value
    
    # Convert angular size to approximate size in plot units
    major_size_rad = np.radians(dso['major_axis_deg'])
    minor_size_rad = np.radians(dso['minor_axis_deg'])
    
    # Create ellipse by plotting points around perimeter
    n_points = 256
    t = np.linspace(0, 2*np.pi, n_points)
    
    # Ellipse in local coordinates
    a = major_size_rad / 2
    b = minor_size_rad / 2
    
    x_local = a * np.cos(t)
    y_local = b * np.sin(t)
    
    # Rotate by position angle
    pa_rad = np.radians(dso['position_angle_deg'])
    cos_pa = np.cos(pa_rad)
    sin_pa = np.sin(pa_rad)
    x_rot = x_local * cos_pa - y_local * sin_pa
    y_rot = x_local * sin_pa + y_local * cos_pa
    
    # Convert to polar plot coordinates using proper spherical trigonometry
    # x_rot and y_rot are angular offsets in radians in the local tangent plane
    # Center point: (azimuth_rad, elevation_rad)
    # Use spherical trigonometry to compute new (azimuth, elevation) for each point
    
    # Convert center to unit vector on sphere
    center_x = np.cos(elevation_rad) * np.cos(azimuth_rad)
    center_y = np.cos(elevation_rad) * np.sin(azimuth_rad)
    center_z = np.sin(elevation_rad)
    center_vec = np.array([center_x, center_y, center_z])
    
    # Create local tangent frame vectors at center
    # North vector (increasing elevation): derivative w.r.t. elevation
    north_x = -np.sin(elevation_rad) * np.cos(azimuth_rad)
    north_y = -np.sin(elevation_rad) * np.sin(azimuth_rad)
    north_z = np.cos(elevation_rad)
    north_vec = np.array([north_x, north_y, north_z])
    
    # East vector (increasing azimuth): derivative w.r.t. azimuth, scaled by cos(elevation)
    east_x = -np.sin(azimuth_rad)
    east_y = np.cos(azimuth_rad)
    east_z = 0.0
    east_vec = np.array([east_x, east_y, east_z])
    
    # For each ellipse point, compute offset vector in tangent plane
    # x_rot is offset in north direction, y_rot is offset in east direction
    # Use small-angle approximation for tangent plane (valid for small angular sizes)
    # Then project back to unit sphere
    
    # Compute offset vectors for all points at once
    offset_vecs = (x_rot[:, np.newaxis] * north_vec + 
                   y_rot[:, np.newaxis] * east_vec)
    
    # Project to unit sphere: new_vec = normalize(center_vec + offset_vec)
    new_vecs = center_vec + offset_vecs
    new_norms = np.linalg.norm(new_vecs, axis=1, keepdims=True)
    new_norms = np.where(new_norms > 1e-10, new_norms, 1.0)  # Avoid division by zero
    new_vecs_normalized = new_vecs / new_norms
    
    # Convert back to spherical coordinates (azimuth, elevation)
    new_x = new_vecs_normalized[:, 0]
    new_y = new_vecs_normalized[:, 1]
    new_z = new_vecs_normalized[:, 2]
    
    new_azimuth = np.arctan2(new_y, new_x)
    new_azimuth = np.mod(new_azimuth, 2 * np.pi)  # Wrap to [0, 2π)
    new_elevation = np.arcsin(np.clip(new_z, -1.0, 1.0))
    
    # Convert to polar plot coordinates
    if is_north_hemisphere:
        ellipse_radial = 0.5 * np.pi - new_elevation
    else:
        ellipse_radial = 0.5 * np.pi + new_elevation
    
    ellipse_azimuth = new_azimuth
    
    # Filter points within valid hemisphere bounds
    if is_north_hemisphere:
        if clip_to_equator:
            valid = ellipse_radial <= equator_radial
        else:
            valid = ellipse_radial >= 0
    else:
        if clip_to_equator:
            valid = ellipse_radial >= equator_radial
        else:
            valid = ellipse_radial >= 0
    
    if np.any(valid):
        # Choose color based on object type
        if dso['object_type'] == 'nebula':
            color = 'cyan'
        elif dso['object_type'] == 'globular':
            color = 'magenta'
        else:  # cluster
            color = 'lime'
        
        # Plot ellipse outline
        ellipse_azimuth_plot = ellipse_azimuth[valid]
        ellipse_radial_plot = ellipse_radial[valid]
        
        ax.plot(ellipse_azimuth_plot, ellipse_radial_plot, 
                color=color, linewidth=1.5, alpha=0.7, transform=ax.transData)
        
        # Add label at center - offset down for magenta (globular), cyan (nebula), and lime (cluster)
        center_in_hemisphere = (is_north_hemisphere and elevation_rad > 0) or (not is_north_hemisphere and elevation_rad < 0)
        if center_in_hemisphere or not clip_to_equator:
            fontsize = _label_font_size(LABEL_FONT_SIZE_RATIO_MED)
            # Use ax.annotate with offset points for proper label positioning
            # Offset down for magenta (globular clusters), cyan (nebulas), and lime (open clusters)
            if color == 'magenta' or color == 'cyan' or color == 'lime':
                offset_y = -10  # 10 points below
            else:
                offset_y = 0  # No offset
            ann = ax.annotate(dso['name'],
                             xy=(azimuth_rad, radial_center),
                             xytext=(0, offset_y),
                             textcoords='offset points',
                             ha='center', va='top' if offset_y < 0 else 'center',
                             color=color, fontsize=fontsize, weight='bold',
                             transform=ax.transData)
            ann.set_path_effects(_text_stroke_effects())

def generate_galactic_hemispheres(target_star_name, search_radius_pc=15, force_refresh=False, star_limit=None, dump_positions=False, magnitude_limit=None, point_size_min=None, point_size_max=None):
    mag_limit = magnitude_limit if magnitude_limit is not None else VISIBLE_MAG_LIMIT
    point_min = point_size_min if point_size_min is not None else STAR_POINT_MIN_SIZE
    point_max = point_size_max if point_size_max is not None else STAR_POINT_MAX_SIZE
    print(f"--- Processing {target_star_name} ---")
    if dump_positions:
        print("Mode: dump positions to DB (no images)")
    
    # Ensure directories exist
    init_database(CACHE_DB)
    IMAGES_DIR.mkdir(parents=True, exist_ok=True)
    
    # Special case: Sol / Sun — observer at solar system barycenter (origin)
    # Positions relative to Sol = positions relative to Earth (no offset)
    if target_star_name.strip().lower() in ('sol', 'sun'):
        print("Target: Sol (Sun). Using origin (0,0,0) — positions = Earth-centric.")
        target_3d = SkyCoord(CartesianRepresentation(0 * u.pc, 0 * u.pc, 0 * u.pc), frame='icrs')
    else:
        # 1. Get Target Star from cache or Gaia
        print(f"Checking cache for target star: {target_star_name}...")
        target_data = get_target_star_from_cache(CACHE_DB, target_star_name)
        
        if target_data is None:
            print(f"Target star not found in cache. Querying Gaia...")
            target_data = get_target_star_from_gaia(target_star_name)
            
            if target_data is None:
                print(f"Error: Could not find target star '{target_star_name}' (cache and Gaia/Simbad lookup failed).", file=sys.stderr)
                print(f"  Check the name (e.g. 'Arcturus' not 'Arcturis') and try again. Skipping this target.", file=sys.stderr)
                return False
            print(f"Caching target star in database...")
            cache_target_star(CACHE_DB, target_data)
            print(f"  Target star cached successfully.")
        else:
            print(f"  Target star found in cache:")
            print(f"    Source ID: {target_data['source_id']}")
            print(f"    RA: {target_data['ra']:.6f}°, Dec: {target_data['dec']:.6f}°")
            print(f"    Parallax: {target_data['parallax']:.4f} mas")
            print(f"    Magnitude: {target_data['phot_g_mean_mag']:.2f}")
        
        # 2. Calculate target star's 3D position relative to Earth
        # Distance formula: d (pc) = 1000 / parallax (mas) — standard parsec definition.
        override = _stardata_get_known_target_parallax_override(target_star_name)
        if override and "parallax_mas" in override:
            target_data["parallax"] = override["parallax_mas"]
            if "ra_deg" in override:
                target_data["ra"] = override["ra_deg"]
            if "dec_deg" in override:
                target_data["dec"] = override["dec_deg"]
            print(f"  Using catalog parallax override for '{target_star_name}': {target_data['parallax']:.4f} mas")
        target_ra = target_data['ra']
        target_dec = target_data['dec']
        target_parallax = target_data['parallax']
        if target_parallax > 0 and np.isfinite(target_parallax):
            target_dist_pc = 1000.0 / target_parallax
            print(f"  Calculated distance from parallax: {target_dist_pc:.2f} pc (~{target_dist_pc * 3.262:.2f} ly)")
        else:
            target_dist_pc = 20.0
            print(f"  Warning: Invalid parallax ({target_parallax}), using default distance: {target_dist_pc} pc")
        
        target_3d = SkyCoord(ra=target_ra*u.deg, dec=target_dec*u.deg, 
                             distance=target_dist_pc*u.pc, frame='icrs')
        print(f"  Target 3D position: RA={target_ra:.6f}°, Dec={target_dec:.6f}°, Distance={target_dist_pc:.2f} pc")
    
    if dump_positions:
        cleared = clear_star_positions_for_target(CACHE_DB, target_star_name)
        print(f"  Cleared {cleared:,} existing rows for target '{target_star_name}' in star_positions_3d.")

    # 3. Ensure cache is populated (only queries Gaia if needed)
    num_stars = ensure_cache_populated(force_refresh=force_refresh, star_limit=star_limit)
    
    if num_stars == 0:
        raise RuntimeError("No stars found in database. Download failed or database is empty.")
    
    print(f"\nProcessing {num_stars:,} stars from SQLite cache...")
    print("(All operations use cached data, not Gaia API)")
    
    # 4. Load all data from SQLite cache and process
    PROCESS_CHUNK_SIZE = 500000
    all_azimuth_rad = []
    all_elevation_rad = []
    all_z = []
    all_m_plot = []
    all_bp_rp = []
    all_distance_pc = []
    
    total_stars = num_stars
    
    if HAS_TQDM:
        pbar = tqdm(total=total_stars, desc='Processing from cache', unit='stars', unit_scale=True)
    else:
        pbar = None
        print(f"Progress: 0/{total_stars:,} stars (0.0%)")
        sys.stdout.flush()
    
    try:
        offset = 0
        while offset < total_stars:
            chunk_data = load_stars_from_cache(CACHE_DB, limit=PROCESS_CHUNK_SIZE, offset=offset)
            if chunk_data is None or len(chunk_data) == 0:
                break
            chunk_size = len(chunk_data)
            
            result = process_star_chunk(chunk_data, target_3d, dump_positions=dump_positions, magnitude_limit=mag_limit)
            
            if result is not None:
                n = len(result['azimuth_rad'])
                all_azimuth_rad.append(result['azimuth_rad'])
                all_elevation_rad.append(result['elevation_rad'])
                all_z.append(result['z'])
                all_m_plot.append(result['m_new'])
                all_bp_rp.append(result['bp_rp'])
                all_distance_pc.append(result['distance_pc'])
                total_valid = sum(len(x) for x in all_azimuth_rad)
                
                if dump_positions:
                    rows = list(zip(
                        np.asarray(result['source_id']).tolist(),
                        result['x_pc'].tolist(),
                        result['y_pc'].tolist(),
                        result['z_pc'].tolist(),
                        result['azimuth_rad'].tolist(),
                        result['elevation_rad'].tolist(),
                        result['m_new'].tolist(),
                    ))
                    insert_star_positions_batch(CACHE_DB, target_star_name, rows)
                
                if HAS_TQDM:
                    pbar.update(chunk_size)
                    pbar.set_postfix({'valid': f'{total_valid:,}'})
                else:
                    percent = ((offset + chunk_size) / total_stars) * 100
                    print(f"  ✓ Processed {chunk_size:,} stars, {n:,} valid (Total: {total_valid:,}, {percent:.1f}%)")
                    sys.stdout.flush()
            else:
                if HAS_TQDM:
                    pbar.update(chunk_size)
                    total_valid = sum(len(x) for x in all_azimuth_rad)
                    pbar.set_postfix({'valid': f'{total_valid:,}'})
                else:
                    percent = ((offset + chunk_size) / total_stars) * 100
                    total_valid = sum(len(x) for x in all_azimuth_rad)
                    print(f"  ✗ No valid stars in chunk (Total: {total_valid:,}, {percent:.1f}%)")
                    sys.stdout.flush()
            
            offset += chunk_size
            del chunk_data, result
    finally:
        if pbar:
            pbar.close()
    
    if not all_azimuth_rad:
        raise RuntimeError("No valid stars found after processing all chunks.")
    
    total_valid = sum(len(x) for x in all_azimuth_rad)
    print(f"\nTotal valid stars: {total_valid:,}")
    
    azimuth_rad = np.concatenate(all_azimuth_rad)
    elevation_rad = np.concatenate(all_elevation_rad)
    z = np.concatenate(all_z)
    m_plot = np.concatenate(all_m_plot)
    bp_rp = np.concatenate(all_bp_rp)
    distance_pc = np.concatenate(all_distance_pc)
    
    # Draw star field in distance order: farthest first so closer stars render on top
    sort_idx = np.argsort(distance_pc)[::-1]
    azimuth_rad = azimuth_rad[sort_idx]
    elevation_rad = elevation_rad[sort_idx]
    z = z[sort_idx]
    m_plot = m_plot[sort_idx]
    bp_rp = bp_rp[sort_idx]
    del distance_pc, all_distance_pc, sort_idx
    
    print(f"  Azimuth range: {np.degrees(np.min(azimuth_rad)):.1f}° to {np.degrees(np.max(azimuth_rad)):.1f}°")
    print(f"  Elevation range: {np.degrees(np.min(elevation_rad)):.1f}° to {np.degrees(np.max(elevation_rad)):.1f}°")
    print(f"  Magnitude range: {np.min(m_plot):.2f} to {np.max(m_plot):.2f}")
    print(f"  Stars with z > 0: {np.sum(z > 0):,}")
    print(f"  Stars with z < 0: {np.sum(z < 0):,}")
    print(f"  Stars with z = 0: {np.sum(z == 0):,}")
    
    # Convert BP-RP color index to RGB colors with magnitude-based transparency.
    # Stars brighter than VISIBLE_MAG_TRANSPARENCY_LIMIT are fully opaque.
    # Stars at mag_limit are rendered at 10% opacity (very faint).
    alpha = np.ones_like(m_plot, dtype=float)
    bright_mask = m_plot <= VISIBLE_MAG_TRANSPARENCY_LIMIT
    faint_mask = m_plot >= mag_limit
    mid_mask = (~bright_mask) & (~faint_mask)

    alpha[bright_mask] = STAR_POINT_MAX_ALPHA
    alpha[faint_mask] = STAR_POINT_MIN_ALPHA
    if np.any(mid_mask):
        # Linearly interpolate opacity between the two thresholds
        span = mag_limit - VISIBLE_MAG_TRANSPARENCY_LIMIT
        frac = (m_plot[mid_mask] - VISIBLE_MAG_TRANSPARENCY_LIMIT) / span
        alpha[mid_mask] = STAR_POINT_MAX_ALPHA - frac * (STAR_POINT_MAX_ALPHA - STAR_POINT_MIN_ALPHA)
    
    # Ensure alpha is fully materialized (not a view) to avoid issues with lazy evaluation.
    # This fixes a bug where increasing VISIBLE_MAG_LIMIT caused fewer stars to be displayed.
    alpha = np.asarray(alpha, dtype=float).copy()

    star_colors = bp_rp_to_rgb(bp_rp, alpha=alpha)
    
    del all_azimuth_rad, all_elevation_rad, all_z, all_m_plot, all_bp_rp
    
    if dump_positions:
        print(f"\nDumped {total_valid:,} star positions to star_positions_3d for target '{target_star_name}'.")
        print("Use SQLite to inspect: SELECT * FROM star_positions_3d WHERE target_star_name = ? LIMIT 20;")
        return
    
    # 5. Plotting Northern and Southern Hemispheres based on z-component
    # Output filenames use target name with spaces replaced by dashes
    filename_base = target_star_name.replace(" ", "-")
    # Scaling for stars
    #
    # IMPORTANT:
    # The shared helper `calculate_point_size_by_magnitude` is used for bright / labelled
    # objects (fixed stars, Sol, DSOs). For the dense star field, we avoid any logic
    # that can ever return size 0 for in-range magnitudes, because:
    # - m_plot has already been filtered by `m_new < VISIBLE_MAG_LIMIT`
    # - runtime logs showed that, for higher VISIBLE_MAG_LIMIT values, nearly all
    #   mid-range stars were being assigned size 0, effectively removing them.
    #
    # Here we therefore apply a simple, explicit piecewise mapping:
    # - 0 <= mag <= 7.0  : full dynamic range from STAR_POINT_MAX_SIZE down to MIN
    # - mag > 7.0       : clamped to STAR_POINT_MIN_SIZE (always > 0)
    # This guarantees that *no* star that passes the magnitude filter will ever have
    # a zero point size, independent of VISIBLE_MAG_LIMIT.
    mag = np.asarray(m_plot, dtype=float)
    point_sizes = np.full_like(mag, point_min, dtype=float)

    mag_threshold = 7.0
    bright_mask = mag <= mag_threshold
    if np.any(bright_mask):
        bright_mag = np.clip(mag[bright_mask], 0.0, mag_threshold)
        bright_mag_norm = bright_mag / mag_threshold  # 0 .. 1
        size_range = point_max - point_min
        point_sizes[bright_mask] = point_max - bright_mag_norm * size_range

    # Get bright galaxies, stars, and deep sky objects
    galaxies = get_bright_galaxies()
    bright_stars = get_bright_stars()
    dso_objects = get_bright_deep_sky_objects()
    
    # Calculate galaxy coordinates from target star's perspective
    galaxy_data = []
    for galaxy in galaxies:
        az, el, z_gal = calculate_galaxy_coordinates(galaxy, target_3d)
        galaxy_data.append({
            'galaxy': galaxy,
            'azimuth_rad': az,
            'elevation_rad': el,
            'z': z_gal
        })
    
    # Calculate Magellanic Clouds coordinates and apparent sizes from target star's perspective
    magellanic_clouds = get_magellanic_clouds()
    mc_data = []
    for mc in magellanic_clouds:
        az, el, z_mc, apparent_major_deg, apparent_minor_deg = calculate_magellanic_cloud_coordinates(mc, target_3d)
        mc_data.append({
            'mc': mc,
            'azimuth_rad': az,
            'elevation_rad': el,
            'z': z_mc,
            'apparent_major_axis_deg': apparent_major_deg,
            'apparent_minor_axis_deg': apparent_minor_deg
        })
    
    # Calculate bright star coordinates and apparent magnitude from target's perspective
    # Skip the target star itself — we're viewing from it, so it shouldn't appear as a labeled point
    target_norm = target_star_name.strip().lower()
    STAR_NAME_ALIASES = {"alpha centauri": "rigil kentaurus"}
    
    def _is_target_star(star):
        star_name_norm = star["name"].lower()
        if target_norm == star_name_norm:
            return True
        canonical = STAR_NAME_ALIASES.get(target_norm)
        return canonical is not None and canonical == star_name_norm
    
    star_data = []
    for star in bright_stars:
        if _is_target_star(star):
            continue
        az, el, z_star = calculate_star_coordinates(star, target_3d)
        
        # Calculate apparent magnitude from target star's perspective
        parallax_mas = star.get('parallax_mas', None)
        if parallax_mas and parallax_mas > 0:
            star_distance_pc = 1000.0 / parallax_mas
        else:
            star_distance_pc = 100000.0  # Very distant
        
        # Get target distance from Earth
        target_cart = target_3d.cartesian.xyz.value
        target_distance_pc = np.linalg.norm(target_cart)
        
        # Calculate distance from target to star
        if star_distance_pc > 100 * target_distance_pc or target_distance_pc < 1e-6:
            # Very distant star - distance from target ≈ distance from Earth
            d_from_target_pc = star_distance_pc
        else:
            # Calculate actual distance from target to star
            star_icrs = SkyCoord(ra=star['ra_deg']*u.deg, dec=star['dec_deg']*u.deg, 
                                distance=star_distance_pc*u.pc, frame='icrs')
            star_cart = star_icrs.cartesian.xyz.value
            vector_from_target = star_cart - target_cart
            d_from_target_pc = np.linalg.norm(vector_from_target)
        
        # Calculate absolute magnitude from Earth's perspective
        # M = m - 5*log10(d/10pc)
        if star_distance_pc > 1e-6:
            M_absolute = star['apparent_mag'] - 5 * np.log10(star_distance_pc / 10.0)
        else:
            M_absolute = star['apparent_mag']
        
        # Calculate apparent magnitude from target's perspective
        # m = M + 5*log10(d/10pc)
        if d_from_target_pc > 1e-6:
            apparent_mag_from_target = M_absolute + 5 * np.log10(d_from_target_pc / 10.0)
        else:
            apparent_mag_from_target = star['apparent_mag']
        
        # Only include stars that are bright enough (magnitude <= mag_limit, same as regular stars)
        if apparent_mag_from_target <= mag_limit:
            star_data.append({
                'star': star,
                'azimuth_rad': az,
                'elevation_rad': el,
                'z': z_star,
                'apparent_mag_from_target': apparent_mag_from_target
            })
    
    # Calculate DSO coordinates
    dso_data = []
    for dso in dso_objects:
        az, el, z_dso = calculate_dso_coordinates(dso, target_3d)
        dso_data.append({
            'dso': dso,
            'azimuth_rad': az,
            'elevation_rad': el,
            'z': z_dso
        })
    
    # Northern Hemisphere (z > 0)
    # Filter out stars with point_sizes <= 0
    valid_size_mask = point_sizes > 0
    north_mask = (z > 0) & valid_size_mask
    if np.any(north_mask):
        fig_north = plt.figure(figsize=(24, 24), facecolor='#000005')
        ax1 = fig_north.add_subplot(111, projection='polar')
        
        # Map elevation to radial distance (center is North Pole, edge is Equator)
        # For north: radial = 90° - elevation (pole at center, equator at edge)
        radial_north = 0.5 * np.pi - elevation_rad[north_mask]
        # Use star colors (already includes alpha), but ensure we have the right shape
        if star_colors.ndim == 2:
            colors_north = star_colors[north_mask]
        else:
            colors_north = star_colors
        ax1.scatter(azimuth_rad[north_mask], radial_north, s=point_sizes[north_mask], c=colors_north, edgecolors='none')
        
        # Plot galaxies that are in or cross the northern hemisphere
        for gd in galaxy_data:
            galaxy_center_elevation = gd['elevation_rad']
            galaxy_major_rad = np.radians(gd['galaxy']['major_axis_deg'] / 2)
            galaxy_top_elevation = galaxy_center_elevation + galaxy_major_rad
            galaxy_bottom_elevation = galaxy_center_elevation - galaxy_major_rad
            
            # Plot if galaxy extends into northern hemisphere
            if galaxy_top_elevation > 0:
                # Check if galaxy crosses equator
                crosses_equator = galaxy_bottom_elevation < 0
                
                if crosses_equator:
                    # Galaxy crosses equator - plot only the northern portion
                    plot_galaxy_on_hemisphere(ax1, gd['galaxy'], gd['azimuth_rad'], 
                                             galaxy_center_elevation, 
                                             is_north_hemisphere=True, 
                                             clip_to_equator=True)
                else:
                    # Entire galaxy is in north hemisphere
                    plot_galaxy_on_hemisphere(ax1, gd['galaxy'], gd['azimuth_rad'], 
                                             galaxy_center_elevation, 
                                             is_north_hemisphere=True, 
                                             clip_to_equator=False)
        
        # Plot Sol as reference star (if target is not Sol)
        if target_star_name.strip().lower() not in ('sol', 'sun'):
            sol_az, sol_el, sol_z = calculate_sol_coordinates(target_3d)
            if sol_z > 0:  # Sol is in northern hemisphere
                plot_sol_reference(ax1, sol_az, sol_el, is_north_hemisphere=True)
        
        # Plot Sagittarius A* as reference point (if target is not Sagittarius A*)
        sgr_a_az, sgr_a_el, sgr_a_z = calculate_sagittarius_a_coordinates(target_3d)
        if sgr_a_z > 0:  # Sagittarius A* is in northern hemisphere
            plot_sagittarius_a_reference(ax1, sgr_a_az, sgr_a_el, is_north_hemisphere=True)
        
        # Plot Magellanic Clouds in or crossing northern hemisphere
        for mc_d in mc_data:
            mc_center_elevation = mc_d['elevation_rad']
            mc_major_rad = np.radians(mc_d['apparent_major_axis_deg'] / 2)
            mc_top_elevation = mc_center_elevation + mc_major_rad
            mc_bottom_elevation = mc_center_elevation - mc_major_rad
            
            if mc_top_elevation > 0:
                crosses_equator = mc_bottom_elevation < 0
                plot_magellanic_cloud(ax1, mc_d['mc'], mc_d['azimuth_rad'], 
                                     mc_center_elevation,
                                     mc_d['apparent_major_axis_deg'],
                                     mc_d['apparent_minor_axis_deg'],
                                     is_north_hemisphere=True,
                                     clip_to_equator=crosses_equator)
        
        # Plot bright stars in northern hemisphere
        for sd in star_data:
            if sd['z'] > 0:  # Star is in northern hemisphere
                plot_star_label(ax1, sd['star'], sd['azimuth_rad'], 
                               sd['elevation_rad'], is_north_hemisphere=True, 
                               apparent_mag_from_target=sd['apparent_mag_from_target'],
                               magnitude_limit=mag_limit,
                               point_size_min=point_min, point_size_max=point_max)
        
        # Plot DSOs in or crossing northern hemisphere
        for dso_d in dso_data:
            dso = dso_d['dso']
            dso_center_elevation = dso_d['elevation_rad']
            dso_major_rad = np.radians(dso['major_axis_deg'] / 2)
            dso_top_elevation = dso_center_elevation + dso_major_rad
            dso_bottom_elevation = dso_center_elevation - dso_major_rad
            
            if dso_top_elevation > 0:
                crosses_equator = dso_bottom_elevation < 0
                plot_dso_on_hemisphere(ax1, dso, dso_d['azimuth_rad'], 
                                       dso_center_elevation, 
                                       is_north_hemisphere=True, 
                                       clip_to_equator=crosses_equator)
        
        title1 = ax1.set_title(f"North Hemisphere\nFrom {target_star_name}", color='white', pad=20,
                               fontsize=_label_font_size(LABEL_FONT_SIZE_RATIO_LARGE))
        title1.set_path_effects(_text_stroke_effects())
        ax1.set_facecolor('#000005')
        ax1.set_yticklabels([]) # Hide radial labels
        tick_labels1 = ax1.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'],
                                           color='gray', fontsize=_label_font_size(LABEL_FONT_SIZE_RATIO_MED))
        for tl in tick_labels1:
            tl.set_path_effects(_text_stroke_effects())
        ax1.grid(True, color='gray', alpha=0.2)
        plt.tight_layout()
        north_image_path = IMAGES_DIR / f"{filename_base}_north_hemisphere.png"
        plt.savefig(north_image_path, facecolor='#000005', bbox_inches='tight', dpi=150)
        plt.close(fig_north)
        print(f"Saved: {north_image_path} ({np.sum(north_mask):,} stars)")
    else:
        print(f"Warning: No stars in north hemisphere (z > 0)")
    
    # Southern Hemisphere (z < 0)
    # Filter out stars with point_sizes <= 0
    south_mask = (z < 0) & valid_size_mask
    if np.any(south_mask):
        fig_south = plt.figure(figsize=(24, 24), facecolor='#000005')
        ax2 = fig_south.add_subplot(111, projection='polar')
        
        # Map elevation to radial distance (center is South Pole, edge is Equator)
        # For south: radial = 90° + elevation (pole at center, equator at edge)
        # Note: elevation is negative for south, so this gives positive radial values
        radial_south = 0.5 * np.pi + elevation_rad[south_mask]
        # Use star colors (already includes alpha), but ensure we have the right shape
        if star_colors.ndim == 2:
            colors_south = star_colors[south_mask]
        else:
            colors_south = star_colors
        ax2.scatter(azimuth_rad[south_mask], radial_south, s=point_sizes[south_mask], c=colors_south, edgecolors='none')
        
        # Plot galaxies that are in or cross the southern hemisphere
        for gd in galaxy_data:
            galaxy_center_elevation = gd['elevation_rad']
            galaxy_major_rad = np.radians(gd['galaxy']['major_axis_deg'] / 2)
            galaxy_top_elevation = galaxy_center_elevation + galaxy_major_rad
            galaxy_bottom_elevation = galaxy_center_elevation - galaxy_major_rad
            
            # Plot if galaxy extends into southern hemisphere
            if galaxy_bottom_elevation < 0:
                # Check if galaxy crosses equator
                crosses_equator = galaxy_top_elevation > 0
                
                if crosses_equator:
                    # Galaxy crosses equator - plot only the southern portion
                    plot_galaxy_on_hemisphere(ax2, gd['galaxy'], gd['azimuth_rad'], 
                                             galaxy_center_elevation, 
                                             is_north_hemisphere=False, 
                                             clip_to_equator=True)
                else:
                    # Entire galaxy is in south hemisphere
                    plot_galaxy_on_hemisphere(ax2, gd['galaxy'], gd['azimuth_rad'], 
                                             galaxy_center_elevation, 
                                             is_north_hemisphere=False, 
                                             clip_to_equator=False)
        
        # Plot Sol as reference star (if target is not Sol)
        if target_star_name.strip().lower() not in ('sol', 'sun'):
            sol_az, sol_el, sol_z = calculate_sol_coordinates(target_3d)
            if sol_z < 0:  # Sol is in southern hemisphere
                plot_sol_reference(ax2, sol_az, sol_el, is_north_hemisphere=False)
        
        # Plot Sagittarius A* as reference point
        sgr_a_az, sgr_a_el, sgr_a_z = calculate_sagittarius_a_coordinates(target_3d)
        if sgr_a_z < 0:  # Sagittarius A* is in southern hemisphere
            plot_sagittarius_a_reference(ax2, sgr_a_az, sgr_a_el, is_north_hemisphere=False)
        
        # Plot Magellanic Clouds in or crossing southern hemisphere
        for mc_d in mc_data:
            mc_center_elevation = mc_d['elevation_rad']
            mc_major_rad = np.radians(mc_d['apparent_major_axis_deg'] / 2)
            mc_top_elevation = mc_center_elevation + mc_major_rad
            mc_bottom_elevation = mc_center_elevation - mc_major_rad
            
            if mc_bottom_elevation < 0:
                crosses_equator = mc_top_elevation > 0
                plot_magellanic_cloud(ax2, mc_d['mc'], mc_d['azimuth_rad'], 
                                     mc_center_elevation,
                                     mc_d['apparent_major_axis_deg'],
                                     mc_d['apparent_minor_axis_deg'],
                                     is_north_hemisphere=False,
                                     clip_to_equator=crosses_equator)
        
        # Plot bright stars in southern hemisphere
        for sd in star_data:
            if sd['z'] < 0:  # Star is in southern hemisphere
                plot_star_label(ax2, sd['star'], sd['azimuth_rad'], 
                               sd['elevation_rad'], is_north_hemisphere=False, 
                               apparent_mag_from_target=sd['apparent_mag_from_target'],
                               magnitude_limit=mag_limit,
                               point_size_min=point_min, point_size_max=point_max)
        
        # Plot DSOs in or crossing southern hemisphere
        for dso_d in dso_data:
            dso = dso_d['dso']
            dso_center_elevation = dso_d['elevation_rad']
            dso_major_rad = np.radians(dso['major_axis_deg'] / 2)
            dso_top_elevation = dso_center_elevation + dso_major_rad
            dso_bottom_elevation = dso_center_elevation - dso_major_rad
            
            if dso_bottom_elevation < 0:
                crosses_equator = dso_top_elevation > 0
                plot_dso_on_hemisphere(ax2, dso, dso_d['azimuth_rad'], 
                                       dso_center_elevation, 
                                       is_north_hemisphere=False, 
                                       clip_to_equator=crosses_equator)
        
        title2 = ax2.set_title(f"South Hemisphere\nFrom {target_star_name}", color='white', pad=20,
                               fontsize=_label_font_size(LABEL_FONT_SIZE_RATIO_LARGE))
        title2.set_path_effects(_text_stroke_effects())
        ax2.set_facecolor('#000005')
        ax2.set_yticklabels([]) # Hide radial labels
        tick_labels2 = ax2.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'],
                                           color='gray', fontsize=_label_font_size(LABEL_FONT_SIZE_RATIO_MED))
        for tl in tick_labels2:
            tl.set_path_effects(_text_stroke_effects())
        ax2.grid(True, color='gray', alpha=0.2)
        plt.tight_layout()
        south_image_path = IMAGES_DIR / f"{filename_base}_south_hemisphere.png"
        plt.savefig(south_image_path, facecolor='#000005', bbox_inches='tight', dpi=150)
        plt.close(fig_south)
        print(f"Saved: {south_image_path} ({np.sum(south_mask):,} stars)")
    else:
        print(f"Warning: No stars in south hemisphere (z < 0)")

    # East / West hemispheres by azimuth, composited into a single image:
    # East:   0°  <= azimuth < 180°
    # West: 180° <= azimuth < 360°
    # Use a full-sky polar projection (north pole at center, south pole at edge),
    # then crop by azimuth range.
    radial_full = 0.5 * np.pi - elevation_rad

    east_mask = (azimuth_rad >= 0.0) & (azimuth_rad < np.pi) & valid_size_mask
    west_mask = (azimuth_rad >= np.pi) & (azimuth_rad < 2 * np.pi) & valid_size_mask

    if np.any(east_mask) or np.any(west_mask):
        fig_east_west = plt.figure(figsize=(24, 24), facecolor='#000005')
        ax_ew = fig_east_west.add_subplot(111, projection='polar')

        # Plot east stars
        if np.any(east_mask):
            if star_colors.ndim == 2:
                colors_east = star_colors[east_mask]
            else:
                colors_east = star_colors

            ax_ew.scatter(
                azimuth_rad[east_mask],
                radial_full[east_mask],
                s=point_sizes[east_mask],
                c=colors_east,
                edgecolors='none',
            )

        # Plot west stars
        if np.any(west_mask):
            if star_colors.ndim == 2:
                colors_west = star_colors[west_mask]
            else:
                colors_west = star_colors

            ax_ew.scatter(
                azimuth_rad[west_mask],
                radial_full[west_mask],
                s=point_sizes[west_mask],
                c=colors_west,
                edgecolors='none',
            )

        # Plot Sol as reference star (if target is not Sol, full sky view)
        if target_star_name.strip().lower() not in ('sol', 'sun'):
            sol_az, sol_el, sol_z = calculate_sol_coordinates(target_3d)
            if sol_el > -np.pi/2 and sol_el < np.pi/2:  # Within visible range
                plot_sol_reference(ax_ew, sol_az, sol_el, is_north_hemisphere=True)
        
        # Plot Sagittarius A* as reference point (full sky view)
        sgr_a_az, sgr_a_el, sgr_a_z = calculate_sagittarius_a_coordinates(target_3d)
        # Use north hemisphere plotting logic for full sky (radial_full uses north convention)
        if sgr_a_el > -np.pi/2 and sgr_a_el < np.pi/2:  # Within visible range
            plot_sagittarius_a_reference(ax_ew, sgr_a_az, sgr_a_el, is_north_hemisphere=True)
        
        # Plot Magellanic Clouds (full sky view)
        for mc_d in mc_data:
            mc_center_elevation = mc_d['elevation_rad']
            # Plot if visible in full sky view
            if mc_center_elevation > -np.pi/2 and mc_center_elevation < np.pi/2:
                plot_magellanic_cloud(ax_ew, mc_d['mc'], mc_d['azimuth_rad'], 
                                     mc_center_elevation,
                                     mc_d['apparent_major_axis_deg'],
                                     mc_d['apparent_minor_axis_deg'],
                                     is_north_hemisphere=True,
                                     clip_to_equator=False)

        title_ew = ax_ew.set_title(f"East / West Hemispheres\nFrom {target_star_name}", color='white', pad=20,
                                   fontsize=_label_font_size(LABEL_FONT_SIZE_RATIO_LARGE))
        title_ew.set_path_effects(_text_stroke_effects())
        ax_ew.set_facecolor('#000005')
        ax_ew.set_yticklabels([])  # Hide radial labels
        tick_labels_ew = ax_ew.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'],
                                               color='gray', fontsize=_label_font_size(LABEL_FONT_SIZE_RATIO_MED))
        for tl in tick_labels_ew:
            tl.set_path_effects(_text_stroke_effects())
        ax_ew.grid(True, color='gray', alpha=0.2)
        plt.tight_layout()
        ew_image_path = IMAGES_DIR / f"{filename_base}_360_degree.png"
        plt.savefig(ew_image_path, facecolor='#000005', bbox_inches='tight', dpi=150)
        plt.close(fig_east_west)
        print(f"Saved: {ew_image_path} ({np.sum(east_mask | west_mask):,} stars)")
    else:
        print("Warning: No stars in east or west hemisphere (0°–360° azimuth)")
    return True

def main():
    """Main entry point for the script."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Generate galactic hemisphere maps from Gaia DR3 data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  {sys.argv[0]} HD118246                    # Generate maps with ~{DEFAULT_STAR_LIMIT:,} stars (default)
  {sys.argv[0]} "Sol,Aldebaran,Proxima Centauri"  # Multiple stars (comma-separated; quote names with spaces)
  {sys.argv[0]} HD118246 100000             # Generate maps with ~100,000 stars
  {sys.argv[0]} Sol 50000 --magnitude-limit 13  # Include stars up to magnitude 13
  {sys.argv[0]} Sol 50000 --dump-positions  # Dump 3D positions for Sol (no images); verify calc
        """
    )
    parser.add_argument(
        'target_star',
        type=str,
        help='Target star name(s), comma-separated. Quote names with spaces (e.g., "Sol,Aldebaran,Proxima Centauri")'
    )
    parser.add_argument(
        'stars',
        type=int,
        nargs='?',
        default=DEFAULT_STAR_LIMIT,
        help=f'Number of stars to download (default: {DEFAULT_STAR_LIMIT:,})'
    )
    parser.add_argument(
        '--force-refresh',
        action='store_true',
        help='Force re-download of cached data'
    )
    parser.add_argument(
        '--dump-positions',
        action='store_true',
        help='Write 3D positions to star_positions_3d table; skip image generation'
    )
    parser.add_argument(
        '--magnitude-limit',
        type=float,
        default=None,
        metavar='MAG',
        help=f'Faintest magnitude to render (default: {VISIBLE_MAG_LIMIT} from constants)'
    )
    parser.add_argument(
        '--point-size-min',
        type=float,
        default=None,
        metavar='SIZE',
        help=f'Smallest star point size in plot (default: {STAR_POINT_MIN_SIZE})'
    )
    parser.add_argument(
        '--point-size-max',
        type=float,
        default=None,
        metavar='SIZE',
        help=f'Largest star point size in plot (default: {STAR_POINT_MAX_SIZE})'
    )
    args = parser.parse_args()
    
    # Parse comma-separated target stars; strip whitespace and surrounding quotes
    target_stars = []
    for part in args.target_star.split(','):
        name = part.strip().strip('"\'')
        if name:
            target_stars.append(name)
    
    if not target_stars:
        parser.error("At least one target star name is required")
    
    failed = []
    for target_star_name in target_stars:
        ok = generate_galactic_hemispheres(
            target_star_name,
            force_refresh=args.force_refresh,
            star_limit=args.stars,
            dump_positions=args.dump_positions,
            magnitude_limit=args.magnitude_limit,
            point_size_min=args.point_size_min,
            point_size_max=args.point_size_max
        )
        if ok is False:
            failed.append(target_star_name)
    if failed:
        print(f"\nSkipped {len(failed)} target(s) due to lookup failure: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)
    return 0

if __name__ == "__main__":
    sys.exit(main())
