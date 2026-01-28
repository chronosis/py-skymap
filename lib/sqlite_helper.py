import sqlite3
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table


def init_database(db_path: Path):
    """Initialize SQLite database with gaia_source and star_positions_3d tables."""
    # Ensure the directory exists
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Create table with source_id as primary key
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS gaia_source (
            source_id INTEGER PRIMARY KEY,
            ra REAL NOT NULL,
            dec REAL NOT NULL,
            parallax REAL NOT NULL,
            phot_g_mean_mag REAL NOT NULL,
            bp_rp REAL
        )
        """
    )
    # Note: PRIMARY KEY automatically creates a unique index, so no manual index needed

    # Add bp_rp column if it doesn't exist (for existing databases)
    cursor.execute(
        """
        SELECT COUNT(*) FROM pragma_table_info('gaia_source') WHERE name='bp_rp'
        """
    )
    if cursor.fetchone()[0] == 0:
        cursor.execute("ALTER TABLE gaia_source ADD COLUMN bp_rp REAL")

    # Table for 3D star positions (used with --dump-positions)
    cursor.execute(
        """
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
        """
    )
    cursor.execute(
        """
        CREATE INDEX IF NOT EXISTS idx_star_positions_target
        ON star_positions_3d(target_star_name)
        """
    )

    conn.commit()
    return conn


def get_star_count(db_path: Path) -> int:
    """Get the number of stars in the database."""
    if not db_path.exists():
        return 0
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM gaia_source")
    count = cursor.fetchone()[0]
    conn.close()
    return count


def get_target_star_from_cache(db_path: Path, target_star_name: str):
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
        cursor.execute(
            """
            SELECT source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp,
                   SQRT(POWER(ra - ?, 2) + POWER(dec - ?, 2)) AS distance
            FROM gaia_source
            WHERE ABS(ra - ?) < ? AND ABS(dec - ?) < ?
            ORDER BY distance
            LIMIT 1
            """,
            (ra_approx, dec_approx, ra_approx, search_radius_deg, dec_approx, search_radius_deg),
        )

        row = cursor.fetchone()
        conn.close()

        if row:
            distance_deg = row[5]
            print(f"  Found star in cache at distance {distance_deg * 3600:.1f} arcseconds")
            return {
                "source_id": row[0],
                "ra": row[1],
                "dec": row[2],
                "parallax": row[3],
                "phot_g_mean_mag": row[4],
            }
    except Exception as e:  # pragma: no cover - defensive logging
        print(f"  Warning: Could not search cache for target star: {e}")

    return None


def cache_target_star(db_path: Path, target_data) -> bool:
    """Cache the target star in the SQLite database."""
    if target_data is None:
        return False

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    try:
        # Use INSERT OR IGNORE to avoid duplicates
        cursor.execute(
            """
            INSERT OR IGNORE INTO gaia_source (source_id, ra, dec, parallax, phot_g_mean_mag)
            VALUES (?, ?, ?, ?, ?)
            """,
            (
                target_data["source_id"],
                target_data["ra"],
                target_data["dec"],
                target_data["parallax"],
                target_data["phot_g_mean_mag"],
            ),
        )
        conn.commit()
        cached = cursor.rowcount > 0
        conn.close()
        return cached
    except Exception as e:  # pragma: no cover - defensive logging
        print(f"  Error caching target star: {e}")
        conn.close()
        return False


def clear_star_positions_for_target(db_path: Path, target_star_name: str) -> int:
    """Remove all star_positions_3d rows for the given target."""
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute(
        "DELETE FROM star_positions_3d WHERE target_star_name = ?",
        (target_star_name,),
    )
    deleted = cursor.rowcount
    conn.commit()
    conn.close()
    return deleted


def insert_star_positions_batch(db_path: Path, target_star_name: str, rows):
    """Insert a batch of (source_id, x_pc, y_pc, z_pc, azimuth_rad, elevation_rad, magnitude)."""
    if not rows:
        return 0
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    data = [
        (r[0], target_star_name, r[1], r[2], r[3], r[4], r[5], r[6])
        for r in rows
    ]
    cursor.executemany(
        """
        INSERT OR REPLACE INTO star_positions_3d
        (source_id, target_star_name, x_pc, y_pc, z_pc, azimuth_rad, elevation_rad, magnitude)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        data,
    )
    n = cursor.rowcount
    conn.commit()
    conn.close()
    return n


def load_stars_from_cache(db_path: Path, limit=None, offset: int = 0):
    """Load stars from SQLite cache.

    This function ONLY reads from the cache, never queries Gaia.
    """
    if not db_path.exists():
        raise RuntimeError(
            f"Cache database {db_path} does not exist. Run ensure_cache_populated() first."
        )

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    if limit is not None:
        cursor.execute(
            """
            SELECT source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp
            FROM gaia_source
            ORDER BY source_id
            LIMIT ? OFFSET ?
            """,
            (limit, offset),
        )
    else:
        cursor.execute(
            """
            SELECT source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp
            FROM gaia_source
            ORDER BY source_id
            LIMIT ? OFFSET ?
            """,
            (1000000, offset),
        )  # Large limit if none specified

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
    bp_rp_values = [row[5] if row[5] is not None else np.nan for row in rows]

    return Table(
        {
            "source_id": source_ids,
            "ra": ra_values * u.deg,
            "dec": dec_values * u.deg,
            "parallax": parallax_values * u.mas,
            "phot_g_mean_mag": mag_values,
            "bp_rp": bp_rp_values,
        }
    )


def check_sky_coverage_bias(db_path: Path):
    """Check if database has biased sky coverage (e.g., only northern hemisphere).

    Returns (has_bias, message) where has_bias is True if bias detected.
    """
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

