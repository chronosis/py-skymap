import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    # Fallback progress function if tqdm not available
    def tqdm(iterable, **kwargs):
        return iterable

CACHE_DB = Path("gaia_cache/gaia_cache.db")
IMAGES_DIR = Path("images")
CHUNK_SIZE = 500000  # Stars per chunk (reduced to avoid timeouts)
DEFAULT_STAR_LIMIT = 50000  # Default number of stars to download

def init_database(db_path):
    """Initialize SQLite database with gaia_source table."""
    # Ensure the directory exists
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    # Create table with source_id as primary key
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS gaia_source (
            source_id INTEGER PRIMARY KEY,
            ra REAL NOT NULL,
            dec REAL NOT NULL,
            parallax REAL NOT NULL,
            phot_g_mean_mag REAL NOT NULL
        )
    """)
    
    # Create index on source_id for faster lookups (though it's already primary key)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_source_id ON gaia_source(source_id)
    """)
    
    # Table for 3D star positions (used with --dump-positions)
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS star_positions_3d (
            source_id INTEGER NOT NULL,
            target_star_name TEXT NOT NULL,
            x_pc REAL NOT NULL,
            y_pc REAL NOT NULL,
            z_pc REAL NOT NULL,
            azimuth_rad REAL NOT NULL,
            elevation_rad REAL NOT NULL,
            magnitude REAL NOT NULL,
            PRIMARY KEY (source_id, target_star_name)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_star_positions_target ON star_positions_3d(target_star_name)
    """)
    
    conn.commit()
    return conn

def get_star_count(db_path):
    """Get the number of stars in the database."""
    if not db_path.exists():
        return 0
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM gaia_source")
    count = cursor.fetchone()[0]
    conn.close()
    return count

def get_target_star_from_cache(db_path, target_star_name):
    """Check if target star is in cache by searching for matching coordinates.
    
    Uses Simbad to get approximate coordinates, then searches cache for stars
    within 1 arcminute. Returns None if not found in cache.
    """
    if not db_path.exists():
        return None
    
    try:
        # Get approximate coordinates from Simbad
        simbad_coord = SkyCoord.from_name(target_star_name)
        ra_approx = simbad_coord.ra.deg
        dec_approx = simbad_coord.dec.deg
        
        # Search in cache for stars near this position (within 1 arcminute = 0.0167 degrees)
        # Use a simple distance calculation: sqrt((ra_diff)^2 + (dec_diff)^2)
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        
        # Search radius in degrees (1 arcminute)
        search_radius_deg = 0.0167
        
        # Find stars within search radius, ordered by distance
        cursor.execute("""
            SELECT source_id, ra, dec, parallax, phot_g_mean_mag,
                   SQRT(POWER(ra - ?, 2) + POWER(dec - ?, 2)) AS distance
            FROM gaia_source
            WHERE ABS(ra - ?) < ? AND ABS(dec - ?) < ?
            ORDER BY distance
            LIMIT 1
        """, (ra_approx, dec_approx, ra_approx, search_radius_deg, dec_approx, search_radius_deg))
        
        row = cursor.fetchone()
        conn.close()
        
        if row:
            distance_deg = row[5]
            print(f"  Found star in cache at distance {distance_deg * 3600:.1f} arcseconds")
            return {
                'source_id': row[0],
                'ra': row[1],
                'dec': row[2],
                'parallax': row[3],
                'phot_g_mean_mag': row[4]
            }
    except Exception as e:
        print(f"  Warning: Could not search cache for target star: {e}")
    
    return None

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
    """Cache the target star in the SQLite database."""
    if target_data is None:
        return False
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        # Use INSERT OR IGNORE to avoid duplicates
        cursor.execute("""
            INSERT OR IGNORE INTO gaia_source (source_id, ra, dec, parallax, phot_g_mean_mag)
            VALUES (?, ?, ?, ?, ?)
        """, (
            target_data['source_id'],
            target_data['ra'],
            target_data['dec'],
            target_data['parallax'],
            target_data['phot_g_mean_mag']
        ))
        conn.commit()
        cached = cursor.rowcount > 0
        conn.close()
        return cached
    except Exception as e:
        print(f"  Error caching target star: {e}")
        conn.close()
        return False

def clear_star_positions_for_target(db_path, target_star_name):
    """Remove all star_positions_3d rows for the given target."""
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("DELETE FROM star_positions_3d WHERE target_star_name = ?", (target_star_name,))
    deleted = cursor.rowcount
    conn.commit()
    conn.close()
    return deleted

def insert_star_positions_batch(db_path, target_star_name, rows):
    """Insert a batch of (source_id, x_pc, y_pc, z_pc, azimuth_rad, elevation_rad, magnitude)."""
    if not rows:
        return 0
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    data = [
        (r[0], target_star_name, r[1], r[2], r[3], r[4], r[5], r[6])
        for r in rows
    ]
    cursor.executemany("""
        INSERT OR REPLACE INTO star_positions_3d
        (source_id, target_star_name, x_pc, y_pc, z_pc, azimuth_rad, elevation_rad, magnitude)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """, data)
    n = cursor.rowcount
    conn.commit()
    conn.close()
    return n

def load_stars_from_cache(db_path, limit=None, offset=0):
    """Load stars from SQLite cache.
    This function ONLY reads from the cache, never queries Gaia."""
    if not db_path.exists():
        raise RuntimeError(f"Cache database {db_path} does not exist. Run ensure_cache_populated() first.")
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    if limit is not None:
        cursor.execute("""
            SELECT source_id, ra, dec, parallax, phot_g_mean_mag
            FROM gaia_source
            ORDER BY source_id
            LIMIT ? OFFSET ?
        """, (limit, offset))
    else:
        cursor.execute("""
            SELECT source_id, ra, dec, parallax, phot_g_mean_mag
            FROM gaia_source
            ORDER BY source_id
            LIMIT ? OFFSET ?
        """, (1000000, offset))  # Large limit if none specified
    
    rows = cursor.fetchall()
    conn.close()
    
    if not rows:
        return None
    
    # Convert to astropy Table
    source_ids = [row[0] for row in rows]
    ra_values = [row[1] for row in rows]
    dec_values = [row[2] for row in rows]
    parallax_values = [row[3] for row in rows]
    mag_values = [row[4] for row in rows]
    
    return Table({
        'source_id': source_ids,
        'ra': ra_values * u.deg,
        'dec': dec_values * u.deg,
        'parallax': parallax_values * u.mas,
        'phot_g_mean_mag': mag_values
    })

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
                # Need more stars, continue downloading
                print(f"Continuing download to reach {star_limit:,} stars...")
                # Continue to download more data
                return _download_from_gaia(star_limit, existing_count)
    
    # Cache is empty or force_refresh, download from Gaia
    return _download_from_gaia(star_limit, 0)

def _download_from_gaia(star_limit=None, existing_count=0):
    """Download data from Gaia API and store in SQLite cache.
    
    ⚠️ IMPORTANT: This is the ONLY function that queries Gaia API directly.
    All other operations (processing, filtering, plotting) use the SQLite cache.
    
    Args:
        star_limit: Maximum number of stars to download (None = no limit)
        existing_count: Number of stars already in cache (for progress tracking)
    
    Returns:
        Total number of unique stars in cache after download
    """
    conn = init_database(CACHE_DB)
    
    if star_limit is not None:
        print(f"Downloading Gaia DR3 data (limit: {star_limit:,} stars)...")
    else:
        print("Downloading Gaia DR3 data...")
        print("(Skipping count query - will download until no more results)")
    
    chunk_num = 0
    offset = 0
    stars_downloaded = 0
    stars_inserted = existing_count  # Track actually inserted (excluding duplicates), start with existing
    
    # Create progress bar
    if HAS_TQDM:
        if star_limit is not None:
            pbar = tqdm(total=star_limit, initial=existing_count, unit='stars', desc='Downloading', unit_scale=True)
        else:
            pbar = tqdm(initial=existing_count, unit='stars', desc='Downloading', unit_scale=True)
    else:
        pbar = None
        if star_limit is not None:
            percent = (stars_inserted / star_limit) * 100 if star_limit > 0 else 0
            print(f"Progress: {stars_inserted:,}/{star_limit:,} stars in cache ({percent:.1f}%)")
        else:
            print(f"Progress: {stars_inserted:,} stars in cache")
        sys.stdout.flush()
    
    try:
        while True:
            # Check if we've reached the star limit
            if star_limit is not None and stars_inserted >= star_limit:
                print(f"\nReached target of {star_limit:,} stars (inserted: {stars_inserted:,})")
                break
            if HAS_TQDM:
                pbar.set_description(f"Downloading chunk {chunk_num + 1}")
            else:
                print(f"Chunk {chunk_num + 1}: Downloading... (Total so far: {stars_inserted:,} unique stars)")
                sys.stdout.flush()
            
            # Gaia TAP uses TOP (not LIMIT). Paginate via OFFSET.
            # CRITICAL: source_id encodes HEALPix level-8 pixels (spatial regions).
            # Ordering by source_id directly would only download stars from certain
            # sky regions (e.g., mostly northern hemisphere).
            # We order by ra, dec to get uniform sky coverage across all declinations.
            # Filter by phot_rp_mean_mag < 12 to get brighter stars and reduce query size.
            query = f"""
            SELECT TOP {CHUNK_SIZE} source_id, ra, dec, parallax, phot_g_mean_mag
            FROM gaiadr3.gaia_source
            WHERE phot_g_mean_mag IS NOT NULL AND phot_rp_mean_mag < 12
            ORDER BY ra, dec
            OFFSET {offset}
            """
            
            # Retry logic for timeouts
            max_retries = 3
            retry_count = 0
            r = None
            
            while retry_count < max_retries:
                try:
                    job = Gaia.launch_job(query)
                    r = job.get_results()
                    break  # Success, exit retry loop
                except Exception as e:
                    retry_count += 1
                    if "timeout" in str(e).lower() or "408" in str(e) or retry_count >= max_retries:
                        if retry_count >= max_retries:
                            print(f"\n  ✗ Failed after {max_retries} retries. Error: {e}")
                            raise
                        else:
                            wait_time = retry_count * 5  # Exponential backoff
                            print(f"  ⚠ Timeout on chunk {chunk_num + 1}, retrying in {wait_time}s... (attempt {retry_count}/{max_retries})")
                            sys.stdout.flush()
                            time.sleep(wait_time)
                    else:
                        raise  # Re-raise if it's not a timeout
            
            # If we got fewer results than requested, we've reached the end
            if len(r) == 0:
                break
            
            # Check if this chunk would exceed the limit
            if star_limit is not None:
                remaining = star_limit - stars_inserted
                if remaining <= 0:
                    break
                # If this chunk would exceed the limit, truncate it
                if len(r) > remaining:
                    r = r[:remaining]
            
            # Insert into SQLite database (using INSERT OR IGNORE to handle duplicates)
            cursor = conn.cursor()
            chunk_inserted = 0
            
            # Extract data from astropy Table
            source_ids = r['source_id'].data
            ra_values = r['ra'].value if hasattr(r['ra'], 'value') else r['ra']
            dec_values = r['dec'].value if hasattr(r['dec'], 'value') else r['dec']
            parallax_values = r['parallax'].value if hasattr(r['parallax'], 'value') else r['parallax']
            mag_values = r['phot_g_mean_mag'].value if hasattr(r['phot_g_mean_mag'], 'value') else r['phot_g_mean_mag']
            
            # Insert in batch
            data_to_insert = [
                (int(sid), float(ra), float(dec), float(par), float(mag))
                for sid, ra, dec, par, mag in zip(source_ids, ra_values, dec_values, parallax_values, mag_values)
            ]
            
            cursor.executemany("""
                INSERT OR IGNORE INTO gaia_source (source_id, ra, dec, parallax, phot_g_mean_mag)
                VALUES (?, ?, ?, ?, ?)
            """, data_to_insert)
            
            chunk_inserted = cursor.rowcount
            conn.commit()
            
            chunk_stars = len(r)
            stars_downloaded += chunk_stars
            stars_inserted += chunk_inserted
            
            # Update progress bar
            duplicates = chunk_stars - chunk_inserted
            if HAS_TQDM:
                if star_limit is not None:
                    # Don't update beyond total
                    remaining = max(0, star_limit - pbar.n)
                    if remaining > 0:
                        pbar.update(min(chunk_inserted, remaining))
                    pbar.set_postfix({
                        'chunk': chunk_num + 1,
                        'new': f'{chunk_inserted:,}', 
                        'dups': f'{duplicates:,}',
                        'total': f'{stars_inserted:,}'
                    })
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
            
            chunk_num += 1
            offset += CHUNK_SIZE
            
            # Check if we've reached the limit after this chunk
            if star_limit is not None and stars_inserted >= star_limit:
                break
            
            # If we got fewer results than requested, we've reached the end
            if chunk_stars < CHUNK_SIZE:
                break
            
            # Clear memory
            del r, job, cursor
    except KeyboardInterrupt:
        print("\n\nDownload interrupted by user.")
        conn.close()
        if HAS_TQDM and pbar is not None:
            pbar.close()
        raise
    except Exception as e:
        print(f"\n\nError during download: {e}")
        conn.close()
        if HAS_TQDM and pbar is not None:
            pbar.close()
        raise
    finally:
        conn.close()
        if HAS_TQDM and pbar is not None:
            pbar.close()
    
    print(f"\nDownload complete. Processed {chunk_num} chunks ({stars_inserted:,} unique stars inserted, {stars_downloaded:,} total downloaded).")
    return stars_inserted

def process_star_chunk(chunk_data, target_3d, dump_positions=False):
    """Process a chunk of stars and return valid stars for plotting or dumping.
    
    Handles both nearby stars (with valid parallax) and background stars
    (with 0 or invalid parallax) by using a large assumed distance.
    
    If dump_positions=True, magnitude filter (m_new < 6.5) is skipped so all
    valid stars are returned for database output.
    """
    # Extract values and handle units properly
    source_ids = np.asarray(chunk_data['source_id'])
    parallax_values = chunk_data['parallax'].value if hasattr(chunk_data['parallax'], 'value') else chunk_data['parallax']
    ra_values = chunk_data['ra'].value if hasattr(chunk_data['ra'], 'value') else chunk_data['ra']
    dec_values = chunk_data['dec'].value if hasattr(chunk_data['dec'], 'value') else chunk_data['dec']
    mag_values = chunk_data['phot_g_mean_mag'].value if hasattr(chunk_data['phot_g_mean_mag'], 'value') else chunk_data['phot_g_mean_mag']
    
    # Filter for valid RA/Dec and magnitude (must be finite, but no magnitude limit)
    valid_coords_mask = (
        np.isfinite(ra_values) & np.isfinite(dec_values) & 
        np.isfinite(mag_values)  # Only check that magnitude is finite, no upper limit
    )
    if not np.any(valid_coords_mask):
        return None
    
    source_ids = source_ids[valid_coords_mask]
    parallax_values = parallax_values[valid_coords_mask]
    ra_values = ra_values[valid_coords_mask]
    dec_values = dec_values[valid_coords_mask]
    mag_values = mag_values[valid_coords_mask]
    
    # Handle parallax: use actual distance if valid, otherwise use large assumed distance for background stars
    # Background stars (parallax <= 0 or invalid) are very far away - use 100,000 pc as assumed distance
    ASSUMED_BACKGROUND_DISTANCE_PC = 100000.0
    BRIGHT_BACKGROUND_ABS_MAG_THRESHOLD = 8.0  # Absolute magnitude threshold for fixed background objects
    
    valid_parallax_mask = (parallax_values > 0) & np.isfinite(parallax_values)
    d_earth_pc = np.where(
        valid_parallax_mask,
        1000.0 / parallax_values,  # Actual distance from parallax
        ASSUMED_BACKGROUND_DISTANCE_PC  # Assumed distance for background stars
    )
    
    # Filter for reasonable distances (avoid infinite or negative)
    valid_distance_mask = (d_earth_pc > 0) & np.isfinite(d_earth_pc) & (d_earth_pc < 1e6)  # Max 1M pc
    if not np.any(valid_distance_mask):
        return None
    
    source_ids = source_ids[valid_distance_mask]
    d_earth_pc = d_earth_pc[valid_distance_mask]
    ra_values = ra_values[valid_distance_mask]
    dec_values = dec_values[valid_distance_mask]
    mag_values = mag_values[valid_distance_mask]
    has_valid_parallax = valid_parallax_mask[valid_distance_mask]
    
    # Identify extremely bright background objects (parallax 0, abs mag <= 8)
    # These should be treated as fixed background points (same position regardless of target)
    background_mask = ~has_valid_parallax
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
    
    # Vector from target to each star (in ICRS cartesian coordinates)
    # For bright background objects, use their Earth-centric direction directly
    # (they appear at the same projected position regardless of target star)
    vectors_from_target = stars_cart - target_cart  # Shape: (N, 3)
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
    
    # Filter for valid target distances
    valid_target_distance_mask = (d_new_pc > 0) & np.isfinite(d_new_pc) & (d_new_pc < 1e6)
    if not np.any(valid_target_distance_mask):
        return None
    
    source_ids = source_ids[valid_target_distance_mask]
    d_new_pc = d_new_pc[valid_target_distance_mask]
    vectors_from_target = vectors_from_target[valid_target_distance_mask]
    ra_values = ra_values[valid_target_distance_mask]
    dec_values = dec_values[valid_target_distance_mask]
    mag_values = mag_values[valid_target_distance_mask]
    d_earth_pc = d_earth_pc[valid_target_distance_mask]
    has_valid_parallax = has_valid_parallax[valid_target_distance_mask]
    bright_background_mask = bright_background_mask[valid_target_distance_mask]
    
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
    
    # Verify we got valid coordinates
    valid_coord_check = (
        np.isfinite(azimuth_rad) & np.isfinite(elevation_rad) &
        np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    )
    if not np.any(valid_coord_check):
        return None
    
    # Calculate apparent magnitude from target star's perspective
    # For stars with valid parallax: calculate absolute magnitude, then apparent from target
    # For bright background objects: use apparent magnitude directly (fixed background position)
    # For other background stars: use apparent magnitude directly (distance change negligible)
    
    m_new = np.zeros_like(mag_values)
    
    # Bright background objects: use apparent magnitude directly (they appear at same position)
    bright_bg_idx = np.where(bright_background_mask)[0]
    if len(bright_bg_idx) > 0:
        m_new[bright_bg_idx] = mag_values[bright_bg_idx]
    
    # Stars with valid parallax: calculate absolute magnitude, then apparent from target
    valid_parallax_idx = np.where(has_valid_parallax)[0]
    if len(valid_parallax_idx) > 0:
        # Absolute magnitude from Earth's perspective
        M_intrinsic = mag_values[valid_parallax_idx] - 5 * np.log10(d_earth_pc[valid_parallax_idx]) + 5
        # Apparent magnitude from target star's perspective
        m_new[valid_parallax_idx] = M_intrinsic + 5 * np.log10(d_new_pc[valid_parallax_idx]) - 5
    
    # Other background stars (not bright enough to be fixed): use apparent magnitude directly
    # For very distant stars, we can approximate by using the Earth's apparent magnitude
    # since the distance difference is negligible compared to their total distance
    other_background_idx = np.where(~has_valid_parallax & ~bright_background_mask)[0]
    if len(other_background_idx) > 0:
        # Background stars are at assumed distance (100k pc), so apparent magnitude from target
        # is approximately the same as from Earth (both are very far from the star)
        # Use Earth's apparent magnitude directly
        m_new[other_background_idx] = mag_values[other_background_idx]
    
    # Filter by apparent magnitude (visible stars) and valid coordinates, unless dumping
    if dump_positions:
        mask = np.isfinite(m_new) & valid_coord_check
    else:
        mask = (m_new < 6.5) & np.isfinite(m_new) & valid_coord_check
    
    if not np.any(mask):
        return None
    
    # Extract coordinates and apply mask (x,y,z in pc from vectors_from_target)
    x_pc = vectors_from_target[:, 0][mask]
    y_pc = vectors_from_target[:, 1][mask]
    z_pc = vectors_from_target[:, 2][mask]
    azimuth_rad_masked = azimuth_rad[mask]
    elevation_rad_masked = elevation_rad[mask]
    z_masked = z[mask]  # unit z-component for hemisphere determination
    m_plot = m_new[mask]
    source_ids_masked = source_ids[mask]
    
    # Final validation: ensure all values are finite
    valid_final = (
        np.isfinite(azimuth_rad_masked) & np.isfinite(elevation_rad_masked) &
        np.isfinite(z_masked) & np.isfinite(m_plot) &
        np.isfinite(x_pc) & np.isfinite(y_pc) & np.isfinite(z_pc)
    )
    if not np.any(valid_final):
        return None
    
    return {
        'source_id': source_ids_masked[valid_final],
        'x_pc': x_pc[valid_final],
        'y_pc': y_pc[valid_final],
        'z_pc': z_pc[valid_final],
        'azimuth_rad': azimuth_rad_masked[valid_final],
        'elevation_rad': elevation_rad_masked[valid_final],
        'z': z_masked[valid_final],  # unit z: positive = north, negative = south
        'm_new': m_plot[valid_final]
    }

def get_bright_galaxies():
    """Return data for the four brightest distant background galaxies visible to naked eye.
    
    Returns list of dicts with: name, ra_deg, dec_deg, major_axis_deg, minor_axis_deg, 
    position_angle_deg, apparent_mag, distance_pc
    Galaxies are too distant for parallax, so we use distance estimates.
    """
    return [
        {
            'name': 'LMC',
            'ra_deg': 80.8938,  # Large Magellanic Cloud
            'dec_deg': -69.7561,
            'major_axis_deg': 10.75,  # Angular size
            'minor_axis_deg': 9.33,
            'position_angle_deg': 0,  # Approximate
            'apparent_mag': 0.9,
            'distance_pc': 50000.0  # ~163 kly
        },
        {
            'name': 'SMC',
            'ra_deg': 13.1867,  # Small Magellanic Cloud
            'dec_deg': -72.8286,
            'major_axis_deg': 5.20,
            'minor_axis_deg': 3.25,
            'position_angle_deg': 45,  # Approximate
            'apparent_mag': 2.7,
            'distance_pc': 60000.0  # ~196 kly
        },
        {
            'name': 'Andromeda',
            'ra_deg': 10.6847,  # M31
            'dec_deg': 41.2687,
            'major_axis_deg': 3.167,  # ~190 arcmin
            'minor_axis_deg': 1.0,     # ~60 arcmin
            'position_angle_deg': 35,  # Approximate
            'apparent_mag': 3.4,
            'distance_pc': 2540000.0  # ~2.54 Mly
        },
        {
            'name': 'Triangulum',
            'ra_deg': 23.4621,  # M33
            'dec_deg': 30.6602,
            'major_axis_deg': 1.0,     # ~60 arcmin
            'minor_axis_deg': 0.7,     # ~40 arcmin
            'position_angle_deg': 23,  # Approximate
            'apparent_mag': 5.7,
            'distance_pc': 3000000.0  # ~3.0 Mly
        }
    ]

def get_bright_stars():
    """Return data for ~100 brightest and most well-known stars.
    
    Returns list of dicts with: name, ra_deg, dec_deg, apparent_mag, parallax_mas
    Parallax values in milliarcseconds (mas). If None or 0, object is too distant for parallax measurement.
    """
    return [
        # Top 20 brightest stars
        # Parallax values from Hipparcos/Gaia (mas). Distance (pc) = 1000 / parallax_mas
        {'name': 'Sirius', 'ra_deg': 101.2872, 'dec_deg': -16.7161, 'apparent_mag': -1.46, 'parallax_mas': 379.21},  # 2.64 pc
        {'name': 'Canopus', 'ra_deg': 95.9880, 'dec_deg': -52.6957, 'apparent_mag': -0.74, 'parallax_mas': 10.43},  # 95.9 pc
        {'name': 'Rigil Kentaurus', 'ra_deg': 219.9009, 'dec_deg': -60.8356, 'apparent_mag': -0.27, 'parallax_mas': 754.81},  # 1.32 pc (Alpha Centauri)
        {'name': 'Arcturus', 'ra_deg': 213.9153, 'dec_deg': 19.1824, 'apparent_mag': -0.05, 'parallax_mas': 88.83},  # 11.26 pc
        {'name': 'Vega', 'ra_deg': 279.2347, 'dec_deg': 38.7837, 'apparent_mag': 0.03, 'parallax_mas': 130.23},  # 7.68 pc
        {'name': 'Capella', 'ra_deg': 79.1723, 'dec_deg': 45.9980, 'apparent_mag': 0.08, 'parallax_mas': 77.29},  # 12.94 pc
        {'name': 'Rigel', 'ra_deg': 78.6345, 'dec_deg': -8.2016, 'apparent_mag': 0.13, 'parallax_mas': 3.78},  # 264.6 pc
        {'name': 'Procyon', 'ra_deg': 114.8255, 'dec_deg': 5.2249, 'apparent_mag': 0.34, 'parallax_mas': 286.05},  # 3.50 pc
        {'name': 'Betelgeuse', 'ra_deg': 88.7929, 'dec_deg': 7.4071, 'apparent_mag': 0.42, 'parallax_mas': 4.51},  # 221.7 pc
        {'name': 'Achernar', 'ra_deg': 24.4285, 'dec_deg': -57.2368, 'apparent_mag': 0.46, 'parallax_mas': 23.39},  # 42.75 pc
        {'name': 'Hadar', 'ra_deg': 210.9559, 'dec_deg': -60.3730, 'apparent_mag': 0.61, 'parallax_mas': 8.32},  # 120.2 pc
        {'name': 'Altair', 'ra_deg': 297.6958, 'dec_deg': 8.8683, 'apparent_mag': 0.76, 'parallax_mas': 194.44},  # 5.14 pc
        {'name': 'Acrux', 'ra_deg': 186.6496, 'dec_deg': -63.0991, 'apparent_mag': 0.77, 'parallax_mas': 10.13},  # 98.7 pc
        {'name': 'Aldebaran', 'ra_deg': 68.9802, 'dec_deg': 16.5093, 'apparent_mag': 0.85, 'parallax_mas': 48.94},  # 20.43 pc
        {'name': 'Antares', 'ra_deg': 247.3519, 'dec_deg': -26.4320, 'apparent_mag': 0.96, 'parallax_mas': 5.89},  # 169.8 pc
        {'name': 'Spica', 'ra_deg': 201.2983, 'dec_deg': -11.1613, 'apparent_mag': 0.98, 'parallax_mas': 12.44},  # 80.39 pc
        {'name': 'Pollux', 'ra_deg': 116.3289, 'dec_deg': 28.0262, 'apparent_mag': 1.14, 'parallax_mas': 96.54},  # 10.36 pc
        {'name': 'Fomalhaut', 'ra_deg': 344.4127, 'dec_deg': -29.6222, 'apparent_mag': 1.16, 'parallax_mas': 130.08},  # 7.69 pc
        {'name': 'Deneb', 'ra_deg': 310.3579, 'dec_deg': 45.2803, 'apparent_mag': 1.25, 'parallax_mas': 2.31},  # 433.0 pc
        {'name': 'Mimosa', 'ra_deg': 191.9303, 'dec_deg': -59.6888, 'apparent_mag': 1.25, 'parallax_mas': 12.59},  # 79.4 pc
        # Additional well-known stars
        {'name': 'Polaris', 'ra_deg': 37.9546, 'dec_deg': 89.2641, 'apparent_mag': 1.98, 'parallax_mas': 7.54},  # 132.6 pc
        {'name': 'Regulus', 'ra_deg': 152.0929, 'dec_deg': 11.9672, 'apparent_mag': 1.35, 'parallax_mas': 42.09},  # 23.76 pc
        {'name': 'Castor', 'ra_deg': 113.6494, 'dec_deg': 31.8883, 'apparent_mag': 1.58, 'parallax_mas': 66.50},  # 15.04 pc
        {'name': 'Bellatrix', 'ra_deg': 81.2828, 'dec_deg': 6.3497, 'apparent_mag': 1.64, 'parallax_mas': 13.42},  # 74.5 pc
        {'name': 'Elnath', 'ra_deg': 81.5728, 'dec_deg': 28.6075, 'apparent_mag': 1.65, 'parallax_mas': 23.84},  # 41.95 pc
        {'name': 'Miaplacidus', 'ra_deg': 138.2999, 'dec_deg': -69.7172, 'apparent_mag': 1.67, 'parallax_mas': 20.71},  # 48.3 pc
        {'name': 'Alnilam', 'ra_deg': 84.0534, 'dec_deg': -1.2019, 'apparent_mag': 1.69, 'parallax_mas': 1.65},  # 606 pc
        {'name': 'Alnitak', 'ra_deg': 85.1897, 'dec_deg': -1.9426, 'apparent_mag': 1.74, 'parallax_mas': 4.43},  # 225.7 pc
        {'name': 'Mirfak', 'ra_deg': 51.0807, 'dec_deg': 49.8612, 'apparent_mag': 1.79, 'parallax_mas': 6.44},  # 155.3 pc
        {'name': 'Dubhe', 'ra_deg': 165.9320, 'dec_deg': 61.7510, 'apparent_mag': 1.79, 'parallax_mas': 26.54},  # 37.68 pc
        {'name': 'Wezen', 'ra_deg': 107.0979, 'dec_deg': -26.3932, 'apparent_mag': 1.83, 'parallax_mas': 2.43},  # 411.5 pc
        {'name': 'Alkaid', 'ra_deg': 206.8852, 'dec_deg': 49.3133, 'apparent_mag': 1.86, 'parallax_mas': 31.88},  # 31.37 pc
        {'name': 'Sargas', 'ra_deg': 255.9867, 'dec_deg': -42.9978, 'apparent_mag': 1.87, 'parallax_mas': 5.89},  # 169.8 pc
        {'name': 'Avior', 'ra_deg': 125.6285, 'dec_deg': -59.5095, 'apparent_mag': 1.86, 'parallax_mas': 7.51},  # 133.2 pc
        {'name': 'Menkalinan', 'ra_deg': 89.8822, 'dec_deg': 44.9474, 'apparent_mag': 1.90, 'parallax_mas': 40.16},  # 24.90 pc
        {'name': 'Atria', 'ra_deg': 252.1662, 'dec_deg': -69.0277, 'apparent_mag': 1.91, 'parallax_mas': 8.33},  # 120.0 pc
        {'name': 'Alhena', 'ra_deg': 99.4279, 'dec_deg': 16.5403, 'apparent_mag': 1.93, 'parallax_mas': 30.49},  # 32.80 pc
        {'name': 'Peacock', 'ra_deg': 306.4119, 'dec_deg': -56.7351, 'apparent_mag': 1.94, 'parallax_mas': 7.56},  # 132.3 pc
        {'name': 'Alsephina', 'ra_deg': 140.5284, 'dec_deg': -54.7088, 'apparent_mag': 1.96, 'parallax_mas': 8.46},  # 118.2 pc
        {'name': 'Mirzam', 'ra_deg': 95.6749, 'dec_deg': -17.9559, 'apparent_mag': 1.98, 'parallax_mas': 3.29},  # 304.0 pc
        {'name': 'Alphard', 'ra_deg': 141.8968, 'dec_deg': -8.6586, 'apparent_mag': 1.99, 'parallax_mas': 41.35},  # 24.18 pc
        {'name': 'Algieba', 'ra_deg': 154.9926, 'dec_deg': 19.8415, 'apparent_mag': 2.01, 'parallax_mas': 25.96},  # 38.52 pc
        {'name': 'Diphda', 'ra_deg': 10.8974, 'dec_deg': -17.9866, 'apparent_mag': 2.04, 'parallax_mas': 33.62},  # 29.74 pc
        {'name': 'Mizar', 'ra_deg': 200.9814, 'dec_deg': 54.9254, 'apparent_mag': 2.04, 'parallax_mas': 41.73},  # 23.96 pc
        {'name': 'Nunki', 'ra_deg': 283.8164, 'dec_deg': -26.2961, 'apparent_mag': 2.05, 'parallax_mas': 13.87},  # 72.1 pc
        {'name': 'Kaus Australis', 'ra_deg': 276.0430, 'dec_deg': -34.3846, 'apparent_mag': 1.79, 'parallax_mas': 8.34},  # 119.9 pc
        {'name': 'Sadr', 'ra_deg': 305.5571, 'dec_deg': 40.2567, 'apparent_mag': 2.23, 'parallax_mas': 1.78},  # 561.8 pc
        {'name': 'Eltanin', 'ra_deg': 262.6082, 'dec_deg': 51.4889, 'apparent_mag': 2.24, 'parallax_mas': 21.09},  # 47.4 pc
        {'name': 'Kaus Media', 'ra_deg': 274.4067, 'dec_deg': -29.8281, 'apparent_mag': 2.70, 'parallax_mas': 8.59},  # 116.4 pc
        {'name': 'Alpheratz', 'ra_deg': 2.0969, 'dec_deg': 29.0904, 'apparent_mag': 2.07, 'parallax_mas': 33.62},  # 29.74 pc
        {'name': 'Mirach', 'ra_deg': 17.4330, 'dec_deg': 35.6206, 'apparent_mag': 2.07, 'parallax_mas': 16.52},  # 60.5 pc
        {'name': 'Rasalgethi', 'ra_deg': 258.6619, 'dec_deg': 14.3903, 'apparent_mag': 2.78, 'parallax_mas': 8.64},  # 115.7 pc
        {'name': 'Kochab', 'ra_deg': 222.6764, 'dec_deg': 74.1555, 'apparent_mag': 2.08, 'parallax_mas': 17.49},  # 57.2 pc
        {'name': 'Saiph', 'ra_deg': 86.9391, 'dec_deg': -9.6696, 'apparent_mag': 2.07, 'parallax_mas': 5.04},  # 198.4 pc
        {'name': 'Hamal', 'ra_deg': 31.7934, 'dec_deg': 23.4624, 'apparent_mag': 2.01, 'parallax_mas': 49.56},  # 20.18 pc
        {'name': 'Algol', 'ra_deg': 47.0422, 'dec_deg': 40.9556, 'apparent_mag': 2.09, 'parallax_mas': 35.14},  # 28.46 pc
        {'name': 'Dschubba', 'ra_deg': 240.0833, 'dec_deg': -22.6217, 'apparent_mag': 2.29, 'parallax_mas': 5.48},  # 182.5 pc
        {'name': 'Zubeneschamali', 'ra_deg': 229.2517, 'dec_deg': -9.3829, 'apparent_mag': 2.61, 'parallax_mas': 18.77},  # 53.3 pc
        {'name': 'Graffias', 'ra_deg': 244.5804, 'dec_deg': -26.4321, 'apparent_mag': 2.62, 'parallax_mas': 5.89},  # 169.8 pc
        {'name': 'Iota Carinae', 'ra_deg': 139.2725, 'dec_deg': -59.2752, 'apparent_mag': 2.21, 'parallax_mas': 4.20},  # 238.1 pc
        {'name': 'Theta Carinae', 'ra_deg': 160.7392, 'dec_deg': -64.3944, 'apparent_mag': 2.74, 'parallax_mas': 5.55},  # 180.2 pc
        {'name': 'Aspidiske', 'ra_deg': 141.5269, 'dec_deg': -64.3944, 'apparent_mag': 2.21, 'parallax_mas': 5.40},  # 185.2 pc
        # Additional bright stars to reach ~100
        {'name': 'Schedar', 'ra_deg': 10.1268, 'dec_deg': 56.5373, 'apparent_mag': 2.24, 'parallax_mas': 14.29},  # 70.0 pc
        {'name': 'Caph', 'ra_deg': 2.2945, 'dec_deg': 59.1498, 'apparent_mag': 2.28, 'parallax_mas': 60.42},  # 16.55 pc
        {'name': 'Achird', 'ra_deg': 9.8322, 'dec_deg': 54.2844, 'apparent_mag': 3.46, 'parallax_mas': 168.45},  # 5.94 pc
        {'name': 'Almach', 'ra_deg': 30.9748, 'dec_deg': 42.3297, 'apparent_mag': 2.10, 'parallax_mas': 16.64},  # 60.1 pc
        {'name': 'Mira', 'ra_deg': 34.8367, 'dec_deg': -2.9774, 'apparent_mag': 2.0, 'parallax_mas': 10.91},  # 91.7 pc
        {'name': 'Menkar', 'ra_deg': 45.5699, 'dec_deg': 4.0897, 'apparent_mag': 2.54, 'parallax_mas': 15.79},  # 63.3 pc
        {'name': 'Baten Kaitos', 'ra_deg': 13.6605, 'dec_deg': -10.3350, 'apparent_mag': 3.74, 'parallax_mas': 16.23},  # 61.6 pc
        {'name': 'Ankaa', 'ra_deg': 6.5708, 'dec_deg': -42.3058, 'apparent_mag': 2.40, 'parallax_mas': 40.90},  # 24.45 pc
        {'name': 'Markab', 'ra_deg': 346.1902, 'dec_deg': 15.2053, 'apparent_mag': 2.49, 'parallax_mas': 15.07},  # 66.4 pc
        {'name': 'Scheat', 'ra_deg': 345.9436, 'dec_deg': 28.0828, 'apparent_mag': 2.44, 'parallax_mas': 16.64},  # 60.1 pc
        {'name': 'Algenib', 'ra_deg': 3.3089, 'dec_deg': 15.1836, 'apparent_mag': 2.83, 'parallax_mas': 12.10},  # 82.6 pc
        {'name': 'Enif', 'ra_deg': 326.0465, 'dec_deg': 9.8750, 'apparent_mag': 2.38, 'parallax_mas': 4.97},  # 201.2 pc
        {'name': 'Homam', 'ra_deg': 340.7508, 'dec_deg': 10.8314, 'apparent_mag': 3.40, 'parallax_mas': 20.23},  # 49.4 pc
        {'name': 'Matar', 'ra_deg': 330.6800, 'dec_deg': 30.2211, 'apparent_mag': 2.99, 'parallax_mas': 13.33},  # 75.0 pc
        {'name': 'Biham', 'ra_deg': 340.3652, 'dec_deg': 6.1978, 'apparent_mag': 3.51, 'parallax_mas': 20.23},  # 49.4 pc
        {'name': 'Sadalbari', 'ra_deg': 344.4127, 'dec_deg': -29.6222, 'apparent_mag': 3.51, 'parallax_mas': 22.46},  # 44.5 pc
    ]

def get_bright_deep_sky_objects():
    """Return data for magnitude 7 or brighter nebulas, globular clusters, and other non-star objects.
    
    Returns list of dicts with: name, ra_deg, dec_deg, major_axis_deg, minor_axis_deg, 
    position_angle_deg, apparent_mag, object_type, distance_pc
    Distance values in parsecs. Where parallax measurements are available, distances are based on those.
    """
    return [
        # Bright Nebulas
        {'name': 'Orion Nebula', 'ra_deg': 83.8221, 'dec_deg': -5.3911, 'major_axis_deg': 1.5, 'minor_axis_deg': 1.0, 'position_angle_deg': 0, 'apparent_mag': 4.0, 'object_type': 'nebula', 'distance_pc': 414.0},  # VLBA parallax
        {'name': 'Carina Nebula', 'ra_deg': 160.8950, 'dec_deg': -59.6856, 'major_axis_deg': 2.0, 'minor_axis_deg': 2.0, 'position_angle_deg': 0, 'apparent_mag': 1.0, 'object_type': 'nebula', 'distance_pc': 2350.0},  # Gaia EDR3
        {'name': 'Eagle Nebula', 'ra_deg': 274.7000, 'dec_deg': -13.8067, 'major_axis_deg': 0.5, 'minor_axis_deg': 0.5, 'position_angle_deg': 0, 'apparent_mag': 6.0, 'object_type': 'nebula', 'distance_pc': 2000.0},  # ~2 kpc
        {'name': 'Lagoon Nebula', 'ra_deg': 271.0000, 'dec_deg': -24.3833, 'major_axis_deg': 1.33, 'minor_axis_deg': 1.0, 'position_angle_deg': 0, 'apparent_mag': 6.0, 'object_type': 'nebula', 'distance_pc': 4100.0},  # ~4.1 kpc
        {'name': 'Trifid Nebula', 'ra_deg': 270.6700, 'dec_deg': -23.0167, 'major_axis_deg': 0.5, 'minor_axis_deg': 0.5, 'position_angle_deg': 0, 'apparent_mag': 6.3, 'object_type': 'nebula', 'distance_pc': 4100.0},  # ~4.1 kpc (near Lagoon)
        {'name': 'Omega Nebula', 'ra_deg': 275.2000, 'dec_deg': -16.1500, 'major_axis_deg': 0.67, 'minor_axis_deg': 0.5, 'position_angle_deg': 0, 'apparent_mag': 6.0, 'object_type': 'nebula', 'distance_pc': 2000.0},  # ~2 kpc
        {'name': 'Rosette Nebula', 'ra_deg': 97.9500, 'dec_deg': 4.9500, 'major_axis_deg': 1.33, 'minor_axis_deg': 1.0, 'position_angle_deg': 0, 'apparent_mag': 4.8, 'object_type': 'nebula', 'distance_pc': 1600.0},  # ~1.6 kpc
        {'name': 'Horsehead Nebula', 'ra_deg': 85.2500, 'dec_deg': -2.4500, 'major_axis_deg': 0.17, 'minor_axis_deg': 0.17, 'position_angle_deg': 0, 'apparent_mag': 6.8, 'object_type': 'nebula', 'distance_pc': 414.0},  # Same region as Orion
        {'name': 'North America Nebula', 'ra_deg': 314.7000, 'dec_deg': 44.0167, 'major_axis_deg': 2.0, 'minor_axis_deg': 1.5, 'position_angle_deg': 0, 'apparent_mag': 4.0, 'object_type': 'nebula', 'distance_pc': 1800.0},  # ~1.8 kpc
        {'name': 'Veil Nebula', 'ra_deg': 312.7500, 'dec_deg': 30.7833, 'major_axis_deg': 3.0, 'minor_axis_deg': 2.5, 'position_angle_deg': 0, 'apparent_mag': 7.0, 'object_type': 'nebula', 'distance_pc': 1470.0},  # ~1.47 kpc (supernova remnant)
        {'name': 'Dumbbell Nebula', 'ra_deg': 299.9017, 'dec_deg': 22.7211, 'major_axis_deg': 0.25, 'minor_axis_deg': 0.17, 'position_angle_deg': 0, 'apparent_mag': 7.4, 'object_type': 'nebula', 'distance_pc': 1360.0},  # Planetary nebula
        {'name': 'Ring Nebula', 'ra_deg': 283.3961, 'dec_deg': 33.0292, 'major_axis_deg': 0.08, 'minor_axis_deg': 0.08, 'position_angle_deg': 0, 'apparent_mag': 8.8, 'object_type': 'nebula', 'distance_pc': 2300.0},  # Planetary nebula
        {'name': 'Helix Nebula', 'ra_deg': 337.4100, 'dec_deg': -20.8333, 'major_axis_deg': 0.5, 'minor_axis_deg': 0.5, 'position_angle_deg': 0, 'apparent_mag': 7.3, 'object_type': 'nebula', 'distance_pc': 695.0},  # ~695 pc (planetary nebula)
        {'name': 'Crab Nebula', 'ra_deg': 83.6331, 'dec_deg': 22.0144, 'major_axis_deg': 0.17, 'minor_axis_deg': 0.17, 'position_angle_deg': 0, 'apparent_mag': 8.4, 'object_type': 'nebula', 'distance_pc': 2000.0},  # ~2 kpc (supernova remnant)
        # Globular Clusters
        {'name': 'Omega Centauri', 'ra_deg': 201.6967, 'dec_deg': -47.4794, 'major_axis_deg': 0.55, 'minor_axis_deg': 0.55, 'position_angle_deg': 0, 'apparent_mag': 3.7, 'object_type': 'globular', 'distance_pc': 5200.0},  # ~5.2 kpc
        {'name': '47 Tucanae', 'ra_deg': 6.0229, 'dec_deg': -72.0814, 'major_axis_deg': 0.5, 'minor_axis_deg': 0.5, 'position_angle_deg': 0, 'apparent_mag': 4.0, 'object_type': 'globular', 'distance_pc': 4500.0},  # ~4.5 kpc
        {'name': 'M13', 'ra_deg': 250.4233, 'dec_deg': 36.4614, 'major_axis_deg': 0.33, 'minor_axis_deg': 0.33, 'position_angle_deg': 0, 'apparent_mag': 5.8, 'object_type': 'globular', 'distance_pc': 7100.0},  # ~7.1 kpc
        {'name': 'M3', 'ra_deg': 205.5483, 'dec_deg': 28.3772, 'major_axis_deg': 0.25, 'minor_axis_deg': 0.25, 'position_angle_deg': 0, 'apparent_mag': 6.2, 'object_type': 'globular', 'distance_pc': 10200.0},  # ~10.2 kpc
        {'name': 'M5', 'ra_deg': 229.6383, 'dec_deg': 2.0811, 'major_axis_deg': 0.25, 'minor_axis_deg': 0.25, 'position_angle_deg': 0, 'apparent_mag': 5.6, 'object_type': 'globular', 'distance_pc': 7500.0},  # ~7.5 kpc
        {'name': 'M4', 'ra_deg': 245.8967, 'dec_deg': -26.5256, 'major_axis_deg': 0.33, 'minor_axis_deg': 0.33, 'position_angle_deg': 0, 'apparent_mag': 5.6, 'object_type': 'globular', 'distance_pc': 2200.0},  # ~2.2 kpc (closest globular)
        {'name': 'M22', 'ra_deg': 279.1000, 'dec_deg': -23.9047, 'major_axis_deg': 0.33, 'minor_axis_deg': 0.33, 'position_angle_deg': 0, 'apparent_mag': 5.1, 'object_type': 'globular', 'distance_pc': 3200.0},  # ~3.2 kpc
        {'name': 'M15', 'ra_deg': 322.4933, 'dec_deg': 12.1672, 'major_axis_deg': 0.25, 'minor_axis_deg': 0.25, 'position_angle_deg': 0, 'apparent_mag': 6.2, 'object_type': 'globular', 'distance_pc': 10400.0},  # ~10.4 kpc
        {'name': 'M2', 'ra_deg': 323.3625, 'dec_deg': -0.8233, 'major_axis_deg': 0.17, 'minor_axis_deg': 0.17, 'position_angle_deg': 0, 'apparent_mag': 6.5, 'object_type': 'globular', 'distance_pc': 11500.0},  # ~11.5 kpc
        {'name': 'M92', 'ra_deg': 259.2817, 'dec_deg': 43.1358, 'major_axis_deg': 0.17, 'minor_axis_deg': 0.17, 'position_angle_deg': 0, 'apparent_mag': 6.4, 'object_type': 'globular', 'distance_pc': 8300.0},  # ~8.3 kpc
        # Open Clusters
        {'name': 'Pleiades', 'ra_deg': 56.8711, 'dec_deg': 24.1053, 'major_axis_deg': 2.0, 'minor_axis_deg': 2.0, 'position_angle_deg': 0, 'apparent_mag': 1.6, 'object_type': 'cluster', 'distance_pc': 135.0},  # Hipparcos/Gaia parallax
        {'name': 'Hyades', 'ra_deg': 66.7325, 'dec_deg': 15.8700, 'major_axis_deg': 5.5, 'minor_axis_deg': 4.0, 'position_angle_deg': 0, 'apparent_mag': 0.5, 'object_type': 'cluster', 'distance_pc': 46.35},  # Hipparcos parallax
        {'name': 'Beehive Cluster', 'ra_deg': 130.1000, 'dec_deg': 19.6667, 'major_axis_deg': 1.5, 'minor_axis_deg': 1.5, 'position_angle_deg': 0, 'apparent_mag': 3.7, 'object_type': 'cluster', 'distance_pc': 187.0},  # M44, ~187 pc
        {'name': 'Double Cluster', 'ra_deg': 34.7500, 'dec_deg': 57.1500, 'major_axis_deg': 1.0, 'minor_axis_deg': 1.0, 'position_angle_deg': 0, 'apparent_mag': 4.3, 'object_type': 'cluster', 'distance_pc': 2300.0},  # h+Chi Persei, ~2.3 kpc
        {'name': 'M6', 'ra_deg': 265.0833, 'dec_deg': -32.2167, 'major_axis_deg': 1.0, 'minor_axis_deg': 1.0, 'position_angle_deg': 0, 'apparent_mag': 4.2, 'object_type': 'cluster', 'distance_pc': 1600.0},  # ~1.6 kpc
        {'name': 'M7', 'ra_deg': 268.4708, 'dec_deg': -34.7917, 'major_axis_deg': 1.33, 'minor_axis_deg': 1.33, 'position_angle_deg': 0, 'apparent_mag': 3.3, 'object_type': 'cluster', 'distance_pc': 980.0},  # ~980 pc
        {'name': 'M11', 'ra_deg': 282.7750, 'dec_deg': -6.2667, 'major_axis_deg': 0.67, 'minor_axis_deg': 0.67, 'position_angle_deg': 0, 'apparent_mag': 5.8, 'object_type': 'cluster', 'distance_pc': 6100.0},  # ~6.1 kpc
        {'name': 'M44', 'ra_deg': 130.1000, 'dec_deg': 19.6667, 'major_axis_deg': 1.5, 'minor_axis_deg': 1.5, 'position_angle_deg': 0, 'apparent_mag': 3.7, 'object_type': 'cluster', 'distance_pc': 187.0},  # Beehive Cluster, same as above
        {'name': 'M45', 'ra_deg': 56.8711, 'dec_deg': 24.1053, 'major_axis_deg': 2.0, 'minor_axis_deg': 2.0, 'position_angle_deg': 0, 'apparent_mag': 1.6, 'object_type': 'cluster', 'distance_pc': 135.0},  # Pleiades, same as above
    ]

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
    
    # Convert to polar plot coordinates
    # For each point, calculate its (azimuth, radial) position
    # Approximate: for small offsets, radial change ≈ x_rot, azimuth change ≈ y_rot / radial_center
    ellipse_radial = radial_center + x_rot
    ellipse_azimuth = azimuth_rad + y_rot / (radial_center + 1e-6)  # Avoid division by zero
    ellipse_azimuth = np.mod(ellipse_azimuth, 2 * np.pi)  # Wrap to [0, 2π)
    
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
            # Calculate text height in radial units
            fontsize = 9
            text_height_radial = calculate_text_height_radial(ax, fontsize, radial_center)
            
            # Position label below center: move down by half text height plus spacing
            if is_north_hemisphere:
                label_radial = radial_center + text_height_radial * 0.5 + 0.01
            else:
                label_radial = radial_center - text_height_radial * 0.5 - 0.01
            
            ax.text(azimuth_rad, label_radial, galaxy['name'], 
                    color='cyan', fontsize=fontsize, ha='center', va='center', weight='bold',
                    transform=ax.transData)

def calculate_galaxy_coordinates(galaxy, target_3d):
    """Calculate azimuth and elevation of a galaxy from target star's perspective.
    
    Uses distance estimates from galaxy data.
    
    Returns: (azimuth_rad, elevation_rad, z)
    """
    # Use distance from galaxy data if available, otherwise use large assumed distance
    galaxy_distance_pc = galaxy.get('distance_pc', 100000.0)
    
    galaxy_icrs = SkyCoord(ra=galaxy['ra_deg']*u.deg, dec=galaxy['dec_deg']*u.deg, 
                           distance=galaxy_distance_pc*u.pc, frame='icrs')
    
    # Get cartesian coordinates from Earth
    galaxy_cart = galaxy_icrs.cartesian.xyz.value  # Shape: (3,)
    target_cart = target_3d.cartesian.xyz.value  # Shape: (3,)
    
    # Vector from target to galaxy
    vector_from_target = galaxy_cart - target_cart  # Shape: (3,)
    
    # Normalize to get unit direction vector
    vector_norm = np.linalg.norm(vector_from_target)
    if vector_norm > 1e-10:
        unit_vector = vector_from_target / vector_norm
    else:
        unit_vector = vector_from_target
    
    # Convert to azimuth and elevation
    x, y, z = unit_vector
    azimuth_rad = np.arctan2(y, x)
    azimuth_rad = np.mod(azimuth_rad, 2 * np.pi)  # Normalize to [0, 2π)
    elevation_rad = np.arcsin(np.clip(z, -1.0, 1.0))
    
    return azimuth_rad, elevation_rad, z

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
    
    # Get target distance from Earth
    target_cart = target_3d.cartesian.xyz.value  # Shape: (3,)
    target_distance_pc = np.linalg.norm(target_cart)
    
    # For stars much farther than the target, the parallax shift is negligible
    # Use Earth-centric direction directly to avoid numerical precision issues
    # For nearby stars, calculate proper 3D position with vector subtraction
    if star_distance_pc > 100 * target_distance_pc or target_distance_pc < 1e-6:
        # Star is much farther than target (or target is at origin) - use Earth-centric direction
        # This is mathematically correct: for D >> d, direction from target ≈ direction from Earth
        star_icrs = SkyCoord(ra=star['ra_deg']*u.deg, dec=star['dec_deg']*u.deg, 
                            distance=1.0*u.pc, frame='icrs')  # Unit distance for direction only
        star_cart = star_icrs.cartesian.xyz.value  # Shape: (3,) - unit direction vector from Earth
        unit_vector = star_cart  # Already normalized, direction from Earth = direction from target
    else:
        # Star is relatively close - calculate proper 3D position with parallax correction
        star_icrs = SkyCoord(ra=star['ra_deg']*u.deg, dec=star['dec_deg']*u.deg, 
                            distance=star_distance_pc*u.pc, frame='icrs')
        
        # Get cartesian coordinates from Earth
        star_cart = star_icrs.cartesian.xyz.value  # Shape: (3,)
        
        # Vector from target to star
        vector_from_target = star_cart - target_cart  # Shape: (3,)
        
        # Normalize to get unit direction vector
        vector_norm = np.linalg.norm(vector_from_target)
        if vector_norm > 1e-10:
            unit_vector = vector_from_target / vector_norm
        else:
            # Star and target are at same position (shouldn't happen for known stars)
            unit_vector = vector_from_target
    
    # Convert to azimuth and elevation
    x, y, z = unit_vector
    azimuth_rad = np.arctan2(y, x)
    azimuth_rad = np.mod(azimuth_rad, 2 * np.pi)  # Normalize to [0, 2π)
    elevation_rad = np.arcsin(np.clip(z, -1.0, 1.0))
    
    return azimuth_rad, elevation_rad, z

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
    
    # Convert direction vectors to galactic frame
    # First, create SkyCoord objects for these directions
    # We'll use a large distance to get the direction
    test_distance = 10000.0 * u.pc
    target_icrs = target_3d
    
    # For each direction, create a point at test_distance and transform to galactic
    # This gives us the direction in galactic coordinates
    intersects_disk = []
    azimuth_result = []
    elevation_result = []
    z_result = []
    
    for i in range(len(x_dir_flat)):
        # Create a point in the direction from target
        dir_vector_icrs = np.array([x_dir_flat[i], y_dir_flat[i], z_dir_flat[i]])
        
        # Get target position in ICRS cartesian
        target_cart_icrs = target_icrs.cartesian.xyz.value
        
        # Point along direction at test distance
        test_point_icrs_cart = target_cart_icrs + dir_vector_icrs * test_distance.value
        
        # Create SkyCoord for this point
        test_point_icrs = SkyCoord(CartesianRepresentation(
            test_point_icrs_cart[0] * u.pc,
            test_point_icrs_cart[1] * u.pc,
            test_point_icrs_cart[2] * u.pc
        ), frame='icrs')
        
        # Transform to galactic
        test_point_gal = test_point_icrs.transform_to('galactic')
        test_cart_gal = test_point_gal.cartesian.xyz.value
        
        # Direction vector in galactic coordinates
        dir_gal = test_cart_gal - target_cart_gal
        dir_gal_norm = np.linalg.norm(dir_gal)
        if dir_gal_norm > 1e-10:
            dir_gal_unit = dir_gal / dir_gal_norm
        else:
            continue
        
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
                  s=50, color='#333333', alpha=0.3, 
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
                           color='#333333', alpha=0.4, edgecolor='#444444', linewidth=0.5,
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
    
    dso_icrs = SkyCoord(ra=dso['ra_deg']*u.deg, dec=dso['dec_deg']*u.deg, 
                        distance=dso_distance_pc*u.pc, frame='icrs')
    
    # Get cartesian coordinates from Earth
    dso_cart = dso_icrs.cartesian.xyz.value  # Shape: (3,)
    target_cart = target_3d.cartesian.xyz.value  # Shape: (3,)
    
    # Vector from target to DSO
    vector_from_target = dso_cart - target_cart  # Shape: (3,)
    
    # Normalize to get unit direction vector
    vector_norm = np.linalg.norm(vector_from_target)
    if vector_norm > 1e-10:
        unit_vector = vector_from_target / vector_norm
    else:
        unit_vector = vector_from_target
    
    # Convert to azimuth and elevation
    x, y, z = unit_vector
    azimuth_rad = np.arctan2(y, x)
    azimuth_rad = np.mod(azimuth_rad, 2 * np.pi)  # Normalize to [0, 2π)
    elevation_rad = np.arcsin(np.clip(z, -1.0, 1.0))
    
    return azimuth_rad, elevation_rad, z

def calculate_text_height_radial(ax, fontsize, radial_position):
    """Calculate the height of text in radial plot coordinates.
    
    Args:
        ax: matplotlib polar axes
        fontsize: Font size in points
        radial_position: Current radial position in plot coordinates
    
    Returns:
        Height in radial units
    """
    # Estimate text height: fontsize in points, convert to pixels
    # 1 point = 4/3 pixels (at 72 DPI), but we'll use a simpler approximation
    # Text height is roughly fontsize pixels
    text_height_pixels = fontsize * 1.2  # Add some padding
    
    # Get the figure and axes dimensions
    fig = ax.figure
    bbox = ax.get_window_extent()
    width_pixels = bbox.width
    height_pixels = bbox.height
    
    # The plot radius in data coordinates is π/2
    # The plot radius in pixels is min(width, height) / 2
    plot_radius_pixels = min(width_pixels, height_pixels) / 2
    
    # Convert pixel height to radial units
    # At the edge (radial = π/2), 1 pixel = (π/2) / plot_radius_pixels
    # But we need to account for the current radial position
    # The scale varies with radial position, but for small offsets we can approximate
    pixel_to_radial = (0.5 * np.pi) / plot_radius_pixels
    
    text_height_radial = text_height_pixels * pixel_to_radial
    
    return text_height_radial

def calculate_sol_coordinates(target_3d):
    """Calculate azimuth and elevation of Sol (Sun) from target star's perspective.
    
    Sol is at the origin (0,0,0) in ICRS coordinates, so the vector from target to Sol
    is simply -target_3d.
    
    Returns: (azimuth_rad, elevation_rad, z)
    """
    # Sol is at origin, so vector from target to Sol is -target_3d
    target_cart = target_3d.cartesian.xyz.value  # Shape: (3,)
    vector_to_sol = -target_cart  # Vector from target to Sol
    
    # Normalize to get unit direction vector
    vector_norm = np.linalg.norm(vector_to_sol)
    if vector_norm > 1e-10:
        unit_vector = vector_to_sol / vector_norm
    else:
        # Target is at origin (Sol itself), return invalid coordinates
        return 0.0, 0.0, 0.0
    
    # Convert to azimuth and elevation
    x, y, z = unit_vector
    azimuth_rad = np.arctan2(y, x)
    azimuth_rad = np.mod(azimuth_rad, 2 * np.pi)  # Normalize to [0, 2π)
    elevation_rad = np.arcsin(np.clip(z, -1.0, 1.0))
    
    return azimuth_rad, elevation_rad, z

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
        
        # Add label
        fontsize = 9
        text_height_radial = calculate_text_height_radial(ax, fontsize, radial)
        
        # Position label fully below the point
        if is_north_hemisphere:
            label_radial = radial + text_height_radial + 0.01
        else:
            label_radial = radial - text_height_radial - 0.01
        
        ax.text(azimuth_rad, label_radial, 'Sol', 
                color='orange', fontsize=fontsize, ha='center', va='top', 
                weight='bold', transform=ax.transData, zorder=10)

def plot_star_label(ax, star, azimuth_rad, elevation_rad, is_north_hemisphere, apparent_mag_from_target):
    """Plot a bright star with label on a polar plot hemisphere.
    
    Args:
        ax: matplotlib polar axes
        star: dict with star data (from get_bright_stars())
        azimuth_rad: azimuth of star in radians
        elevation_rad: elevation of star in radians
        is_north_hemisphere: True for north, False for south
        apparent_mag_from_target: Apparent magnitude from target star's perspective
    """
    # Calculate radial position on polar plot
    if is_north_hemisphere:
        radial = 0.5 * np.pi - elevation_rad  # Pole at center
    else:
        radial = 0.5 * np.pi + elevation_rad  # Pole at center, elevation is negative
    
    # Only plot if in this hemisphere
    if (is_north_hemisphere and elevation_rad > 0) or (not is_north_hemisphere and elevation_rad < 0):
        # Calculate point size based on apparent magnitude from target
        point_size = (7 - apparent_mag_from_target) ** 2.5
        point_size = max(10, min(200, point_size))  # Clamp between 10 and 200
        
        # Plot star as a point
        ax.scatter(azimuth_rad, radial, s=point_size, color='yellow', 
                  alpha=0.9, edgecolors='orange', linewidths=1, transform=ax.transData)
        
        # Calculate text height in radial units
        fontsize = 7
        text_height_radial = calculate_text_height_radial(ax, fontsize, radial)
        
        # Position label fully below the star point
        # For north hemisphere: increase radial (move toward equator/edge)
        # For south hemisphere: decrease radial (move toward pole/center)
        if is_north_hemisphere:
            # Move down: increase radial by text height plus some spacing
            label_radial = radial + text_height_radial + 0.01
        else:
            # Move down: decrease radial by text height plus some spacing
            label_radial = radial - text_height_radial - 0.01
        
        ax.text(azimuth_rad, label_radial, star['name'], 
                color='yellow', fontsize=fontsize, ha='center', va='top', weight='bold',
                transform=ax.transData)

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
    
    # Convert to polar plot coordinates
    ellipse_radial = radial_center + x_rot
    ellipse_azimuth = azimuth_rad + y_rot / (radial_center + 1e-6)
    ellipse_azimuth = np.mod(ellipse_azimuth, 2 * np.pi)
    
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
            # Calculate text height in radial units
            fontsize = 8
            text_height_radial = calculate_text_height_radial(ax, fontsize, radial_center)
            
            # Offset down for magenta (globular clusters), cyan (nebulas), and lime (open clusters)
            if color == 'magenta' or color == 'cyan' or color == 'lime':
                # Position label below center: move down by half text height plus spacing
                if is_north_hemisphere:
                    label_radial = radial_center + text_height_radial * 0.5 + 0.01
                else:
                    label_radial = radial_center - text_height_radial * 0.5 - 0.01
            else:
                label_radial = radial_center
            ax.text(azimuth_rad, label_radial, dso['name'], 
                    color=color, fontsize=fontsize, ha='center', va='center', weight='bold',
                    transform=ax.transData)

def generate_galactic_hemispheres(target_star_name, search_radius_pc=15, force_refresh=False, star_limit=None, dump_positions=False):
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
                raise RuntimeError(f"Could not find target star '{target_star_name}' in Gaia database.")
            
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
        target_ra = target_data['ra']
        target_dec = target_data['dec']
        target_parallax = target_data['parallax']
        
        if target_parallax > 0 and np.isfinite(target_parallax):
            target_dist_pc = 1000.0 / target_parallax
            print(f"  Calculated distance from parallax: {target_dist_pc:.2f} pc")
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
            
            result = process_star_chunk(chunk_data, target_3d, dump_positions=dump_positions)
            
            if result is not None:
                n = len(result['azimuth_rad'])
                all_azimuth_rad.append(result['azimuth_rad'])
                all_elevation_rad.append(result['elevation_rad'])
                all_z.append(result['z'])
                all_m_plot.append(result['m_new'])
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
    
    print(f"  Azimuth range: {np.degrees(np.min(azimuth_rad)):.1f}° to {np.degrees(np.max(azimuth_rad)):.1f}°")
    print(f"  Elevation range: {np.degrees(np.min(elevation_rad)):.1f}° to {np.degrees(np.max(elevation_rad)):.1f}°")
    print(f"  Magnitude range: {np.min(m_plot):.2f} to {np.max(m_plot):.2f}")
    print(f"  Stars with z > 0: {np.sum(z > 0):,}")
    print(f"  Stars with z < 0: {np.sum(z < 0):,}")
    print(f"  Stars with z = 0: {np.sum(z == 0):,}")
    
    del all_azimuth_rad, all_elevation_rad, all_z, all_m_plot
    
    if dump_positions:
        print(f"\nDumped {total_valid:,} star positions to star_positions_3d for target '{target_star_name}'.")
        print("Use SQLite to inspect: SELECT * FROM star_positions_3d WHERE target_star_name = ? LIMIT 20;")
        return
    
    # 5. Plotting Northern and Southern Hemispheres based on z-component
    # Scaling for stars
    point_sizes = (7 - m_plot) ** 2.5

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
    
    # Calculate bright star coordinates and apparent magnitude from target's perspective
    star_data = []
    for star in bright_stars:
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
        
        # Only include stars that are bright enough (magnitude <= 6.5, same as regular stars)
        if apparent_mag_from_target <= 6.5:
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
    north_mask = z > 0
    if np.any(north_mask):
        fig_north = plt.figure(figsize=(24, 24), facecolor='#000005')
        ax1 = fig_north.add_subplot(111, projection='polar')
        
        # Map elevation to radial distance (center is North Pole, edge is Equator)
        # For north: radial = 90° - elevation (pole at center, equator at edge)
        radial_north = 0.5 * np.pi - elevation_rad[north_mask]
        ax1.scatter(azimuth_rad[north_mask], radial_north, s=point_sizes[north_mask], color='white', alpha=0.8)
        
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
        
        # Plot bright stars in northern hemisphere
        for sd in star_data:
            if sd['z'] > 0:  # Star is in northern hemisphere
                plot_star_label(ax1, sd['star'], sd['azimuth_rad'], 
                               sd['elevation_rad'], is_north_hemisphere=True, 
                               apparent_mag_from_target=sd['apparent_mag_from_target'])
        
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
        
        ax1.set_title(f"North Hemisphere\nFrom {target_star_name}", color='white', pad=20)
        ax1.set_facecolor('#000005')
        ax1.set_yticklabels([]) # Hide radial labels
        ax1.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'], color='gray', fontsize=8)
        ax1.grid(True, color='gray', alpha=0.2)
        plt.tight_layout()
        north_image_path = IMAGES_DIR / f"{target_star_name}_north_hemisphere.png"
        plt.savefig(north_image_path, facecolor='#000005', bbox_inches='tight', dpi=150)
        plt.close(fig_north)
        print(f"Saved: {north_image_path} ({np.sum(north_mask):,} stars)")
    else:
        print(f"Warning: No stars in north hemisphere (z > 0)")
    
    # Southern Hemisphere (z < 0)
    south_mask = z < 0
    if np.any(south_mask):
        fig_south = plt.figure(figsize=(24, 24), facecolor='#000005')
        ax2 = fig_south.add_subplot(111, projection='polar')
        
        # Map elevation to radial distance (center is South Pole, edge is Equator)
        # For south: radial = 90° + elevation (pole at center, equator at edge)
        # Note: elevation is negative for south, so this gives positive radial values
        radial_south = 0.5 * np.pi + elevation_rad[south_mask]
        ax2.scatter(azimuth_rad[south_mask], radial_south, s=point_sizes[south_mask], color='white', alpha=0.8)
        
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
        
        # Plot bright stars in southern hemisphere
        for sd in star_data:
            if sd['z'] < 0:  # Star is in southern hemisphere
                plot_star_label(ax2, sd['star'], sd['azimuth_rad'], 
                               sd['elevation_rad'], is_north_hemisphere=False, 
                               apparent_mag_from_target=sd['apparent_mag_from_target'])
        
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
        
        ax2.set_title(f"South Hemisphere\nFrom {target_star_name}", color='white', pad=20)
        ax2.set_facecolor('#000005')
        ax2.set_yticklabels([]) # Hide radial labels
        ax2.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'], color='gray', fontsize=8)
        ax2.grid(True, color='gray', alpha=0.2)
        plt.tight_layout()
        south_image_path = IMAGES_DIR / f"{target_star_name}_south_hemisphere.png"
        plt.savefig(south_image_path, facecolor='#000005', bbox_inches='tight', dpi=150)
        plt.close(fig_south)
        print(f"Saved: {south_image_path} ({np.sum(south_mask):,} stars)")
    else:
        print(f"Warning: No stars in south hemisphere (z < 0)")

def main():
    """Main entry point for the script."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Generate galactic hemisphere maps from Gaia DR3 data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  {sys.argv[0]} HD118246                    # Generate maps with ~{DEFAULT_STAR_LIMIT:,} stars (default)
  {sys.argv[0]} HD118246 100000             # Generate maps with ~100,000 stars
  {sys.argv[0]} Sol 50000 --dump-positions  # Dump 3D positions for Sol (no images); verify calc
        """
    )
    parser.add_argument(
        'target_star',
        type=str,
        help='Target star name (e.g., HD118246)'
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
    args = parser.parse_args()
    
    generate_galactic_hemispheres(
        args.target_star, 
        force_refresh=args.force_refresh,
        star_limit=args.stars,
        dump_positions=args.dump_positions
    )

if __name__ == "__main__":
    main()
