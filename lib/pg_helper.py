import math
import datetime

import numpy as np
import psycopg2
import psycopg2.extras
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table


def _get_conn(dsn: str):
    """Open a new PostgreSQL connection."""
    return psycopg2.connect(dsn)


def _column_exists(cursor, table: str, column: str) -> bool:
    """Check whether *column* exists on *table* via information_schema."""
    cursor.execute(
        """
        SELECT 1 FROM information_schema.columns
        WHERE table_name = %s AND column_name = %s
        LIMIT 1
        """,
        (table, column),
    )
    return cursor.fetchone() is not None


def init_database(dsn: str):
    """Create tables if they don't exist and return an open connection.

    Mirrors the schema of sqlite_helper.init_database but targets PostgreSQL.
    Uses BIGINT for source_id (Gaia source IDs can exceed 32-bit range) and
    DOUBLE PRECISION for coordinate / magnitude columns.
    """
    conn = _get_conn(dsn)
    cursor = conn.cursor()

    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS gaia_source (
            source_id BIGINT PRIMARY KEY,
            ra DOUBLE PRECISION NOT NULL,
            dec DOUBLE PRECISION NOT NULL,
            parallax DOUBLE PRECISION NOT NULL,
            phot_g_mean_mag DOUBLE PRECISION NOT NULL,
            bp_rp DOUBLE PRECISION,
            source TEXT NOT NULL DEFAULT 'gaia'
        )
        """
    )

    if not _column_exists(cursor, "gaia_source", "bp_rp"):
        cursor.execute("ALTER TABLE gaia_source ADD COLUMN bp_rp DOUBLE PRECISION")

    if not _column_exists(cursor, "gaia_source", "source"):
        cursor.execute(
            "ALTER TABLE gaia_source ADD COLUMN source TEXT NOT NULL DEFAULT 'gaia'"
        )
        cursor.execute("UPDATE gaia_source SET source = 'gaia' WHERE source IS NULL")

    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS star_cartesian_earth (
            source_id BIGINT PRIMARY KEY,
            x_pc DOUBLE PRECISION NOT NULL,
            y_pc DOUBLE PRECISION NOT NULL,
            z_pc DOUBLE PRECISION NOT NULL,
            source TEXT NOT NULL DEFAULT 'gaia'
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS simbad_cache (
            main_id TEXT PRIMARY KEY,
            ra DOUBLE PRECISION NOT NULL,
            dec DOUBLE PRECISION NOT NULL,
            parallax_mas DOUBLE PRECISION,
            vmag DOUBLE PRECISION,
            otype TEXT,
            cached_at TEXT NOT NULL
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS simbad_cache_aliases (
            alias TEXT PRIMARY KEY,
            main_id TEXT NOT NULL,
            FOREIGN KEY (main_id) REFERENCES simbad_cache(main_id)
        )
        """
    )

    conn.commit()
    return conn


def get_star_count(dsn: str) -> int:
    """Total number of stars across all sources."""
    conn = _get_conn(dsn)
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM gaia_source")
    count = cursor.fetchone()[0]
    conn.close()
    return count


def get_star_count_by_source(dsn: str, source: str) -> int:
    """Number of stars for a given source ('gaia', 'vizier', or 'simbad')."""
    conn = _get_conn(dsn)
    cursor = conn.cursor()
    cursor.execute(
        "SELECT COUNT(*) FROM gaia_source WHERE source = %s",
        (source,),
    )
    count = cursor.fetchone()[0]
    conn.close()
    return count


def get_target_star_from_cache(dsn: str, target_star_name: str):
    """Search the cache for a star near the Simbad-resolved position.

    Returns a dict with source_id / ra / dec / parallax / phot_g_mean_mag,
    or None if nothing is within 1 arcminute.
    """
    try:
        simbad_coord = SkyCoord.from_name(target_star_name)
        ra_approx = simbad_coord.ra.deg
        dec_approx = simbad_coord.dec.deg

        conn = _get_conn(dsn)
        cursor = conn.cursor()

        search_radius_deg = 0.0167

        cursor.execute(
            """
            SELECT source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp,
                   SQRT(POWER(ra - %s, 2) + POWER(dec - %s, 2)) AS distance
            FROM gaia_source
            WHERE ABS(ra - %s) < %s AND ABS(dec - %s) < %s
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
    except Exception as e:
        print(f"  Warning: Could not search cache for target star: {e}")

    return None


def cache_target_star(dsn: str, target_data) -> bool:
    """Insert the target star into the cache (skip if already present)."""
    if target_data is None:
        return False

    conn = _get_conn(dsn)
    cursor = conn.cursor()

    try:
        cursor.execute(
            """
            INSERT INTO gaia_source
            (source_id, ra, dec, parallax, phot_g_mean_mag, source)
            VALUES (%s, %s, %s, %s, %s, 'gaia')
            ON CONFLICT (source_id) DO NOTHING
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
    except Exception as e:
        print(f"  Error caching target star: {e}")
        conn.close()
        return False


def ra_dec_distance_to_cartesian(
    ra_deg: float, dec_deg: float, d_pc: float
) -> tuple[float, float, float]:
    """Convert ra (deg), dec (deg), distance (pc) to ICRS Cartesian x, y, z in parsecs."""
    ra_rad = math.radians(ra_deg)
    dec_rad = math.radians(dec_deg)
    cos_dec = math.cos(dec_rad)
    x = d_pc * cos_dec * math.cos(ra_rad)
    y = d_pc * cos_dec * math.sin(ra_rad)
    z = d_pc * math.sin(dec_rad)
    return (x, y, z)


def insert_star_cartesian_earth_batch(
    dsn: str,
    rows: list[tuple[int, float, float, float, str]],
) -> int:
    """Bulk-insert into star_cartesian_earth using execute_values for speed."""
    if not rows:
        return 0
    conn = _get_conn(dsn)
    cursor = conn.cursor()
    psycopg2.extras.execute_values(
        cursor,
        """
        INSERT INTO star_cartesian_earth (source_id, x_pc, y_pc, z_pc, source)
        VALUES %s
        ON CONFLICT (source_id) DO UPDATE SET
            x_pc = EXCLUDED.x_pc,
            y_pc = EXCLUDED.y_pc,
            z_pc = EXCLUDED.z_pc,
            source = EXCLUDED.source
        """,
        rows,
        page_size=5000,
    )
    n = cursor.rowcount
    conn.commit()
    conn.close()
    return n


def load_stars_from_cache(dsn: str, limit=None, offset: int = 0):
    """Load stars from PostgreSQL cache (read-only, never queries Gaia).

    Returns an astropy Table with the same columns as the SQLite version.
    """
    conn = _get_conn(dsn)
    cursor = conn.cursor()

    limit_val = limit if limit is not None else 1_000_000
    cursor.execute(
        """
        SELECT g.source_id, g.ra, g.dec, g.parallax, g.phot_g_mean_mag, g.bp_rp,
               c.x_pc, c.y_pc, c.z_pc
        FROM gaia_source g
        LEFT JOIN star_cartesian_earth c ON g.source_id = c.source_id
        ORDER BY g.source_id
        LIMIT %s OFFSET %s
        """,
        (limit_val, offset),
    )

    rows = cursor.fetchall()
    conn.close()

    if not rows:
        return None

    source_ids = [row[0] for row in rows]
    ra_values = [row[1] for row in rows]
    dec_values = [row[2] for row in rows]
    parallax_values = [row[3] for row in rows]
    mag_values = [row[4] for row in rows]
    bp_rp_values = [row[5] if row[5] is not None else np.nan for row in rows]
    x_pc_values = [row[6] if row[6] is not None else np.nan for row in rows]
    y_pc_values = [row[7] if row[7] is not None else np.nan for row in rows]
    z_pc_values = [row[8] if row[8] is not None else np.nan for row in rows]

    return Table(
        {
            "source_id": source_ids,
            "ra": ra_values * u.deg,
            "dec": dec_values * u.deg,
            "parallax": parallax_values * u.mas,
            "phot_g_mean_mag": mag_values,
            "bp_rp": bp_rp_values,
            "x_earth_pc": x_pc_values,
            "y_earth_pc": y_pc_values,
            "z_earth_pc": z_pc_values,
        }
    )


def get_simbad_from_cache(dsn: str, identifier: str):
    """Return a cached Simbad object by identifier (main_id or alias; case-insensitive)."""
    conn = _get_conn(dsn)
    cursor = conn.cursor()
    key = identifier.strip()

    cursor.execute(
        "SELECT main_id FROM simbad_cache_aliases WHERE LOWER(TRIM(alias)) = LOWER(%s) LIMIT 1",
        (key,),
    )
    alias_row = cursor.fetchone()
    if alias_row:
        key = alias_row[0]

    cursor.execute(
        """
        SELECT main_id, ra, dec, parallax_mas, vmag, otype, cached_at
        FROM simbad_cache
        WHERE LOWER(TRIM(main_id)) = LOWER(%s)
        LIMIT 1
        """,
        (key,),
    )
    row = cursor.fetchone()
    conn.close()

    if row is None:
        return None
    return {
        "main_id": row[0],
        "ra": row[1],
        "dec": row[2],
        "parallax_mas": row[3],
        "vmag": row[4],
        "otype": row[5],
        "cached_at": row[6],
    }


def put_simbad_in_cache(
    dsn: str,
    rows: list[dict],
    *,
    requested_identifiers: list[str] | None = None,
) -> int:
    """Insert or update Simbad objects in the cache."""
    if not rows:
        return 0
    now = datetime.datetime.utcnow().isoformat() + "Z"
    conn = _get_conn(dsn)
    cursor = conn.cursor()

    data = [
        (
            r["main_id"],
            float(r["ra"]),
            float(r["dec"]),
            r.get("parallax_mas"),
            r.get("vmag"),
            r.get("otype"),
            now,
        )
        for r in rows
    ]
    psycopg2.extras.execute_values(
        cursor,
        """
        INSERT INTO simbad_cache
        (main_id, ra, dec, parallax_mas, vmag, otype, cached_at)
        VALUES %s
        ON CONFLICT (main_id) DO UPDATE SET
            ra = EXCLUDED.ra,
            dec = EXCLUDED.dec,
            parallax_mas = EXCLUDED.parallax_mas,
            vmag = EXCLUDED.vmag,
            otype = EXCLUDED.otype,
            cached_at = EXCLUDED.cached_at
        """,
        data,
        page_size=1000,
    )
    n = cursor.rowcount

    if requested_identifiers and len(requested_identifiers) >= len(rows):
        for i, r in enumerate(rows):
            alias = requested_identifiers[i].strip()
            if alias and alias.lower() != r["main_id"].lower():
                cursor.execute(
                    """
                    INSERT INTO simbad_cache_aliases (alias, main_id)
                    VALUES (%s, %s)
                    ON CONFLICT (alias) DO UPDATE SET main_id = EXCLUDED.main_id
                    """,
                    (alias, r["main_id"]),
                )

    conn.commit()
    conn.close()
    return n


def get_all_simbad_cached(dsn: str):
    """Return all cached Simbad objects as a list of dicts."""
    conn = _get_conn(dsn)
    cursor = conn.cursor()
    cursor.execute(
        """
        SELECT main_id, ra, dec, parallax_mas, vmag, otype, cached_at
        FROM simbad_cache
        ORDER BY main_id
        """
    )
    rows = cursor.fetchall()
    conn.close()
    return [
        {
            "main_id": r[0],
            "ra": r[1],
            "dec": r[2],
            "parallax_mas": r[3],
            "vmag": r[4],
            "otype": r[5],
            "cached_at": r[6],
        }
        for r in rows
    ]


def check_sky_coverage_bias(dsn: str):
    """Check if database has biased sky coverage (e.g., only northern hemisphere).

    Returns (has_bias, message) where has_bias is True if bias detected.
    """
    conn = _get_conn(dsn)
    cursor = conn.cursor()

    cursor.execute("SELECT COUNT(*) FROM gaia_source WHERE dec < 0")
    n_neg = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM gaia_source")
    total = cursor.fetchone()[0]
    conn.close()

    if total == 0:
        return False, None

    neg_ratio = n_neg / total
    if neg_ratio < 0.1:
        return True, (
            f"WARNING: Database appears to have biased sky coverage. "
            f"Only {n_neg:,} stars ({100*neg_ratio:.1f}%) have dec < 0 (should be ~50%). "
            f"This is likely due to ORDER BY source_id in old downloads. "
            f"Re-download with --force-refresh to get uniform coverage."
        )

    return False, None
