import time
from pathlib import Path

from astroquery.gaia import Gaia
from astropy.table import Table
import numpy as np

from .constants import GAIA_MAX_RETRIES
from .sqlite_helper import init_database
from .progress import HAS_TQDM, tqdm


def _extract_gaia_data_from_table(table: Table):
    """Extract data from Astropy table into list of tuples for database insertion.

    Args:
        table: Astropy Table with columns: source_id, ra, dec, parallax,
               phot_g_mean_mag, bp_rp, phot_bp_mean_mag, phot_rp_mean_mag

    Returns:
        List of tuples: [(source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp), ...]
    """

    def get_value(col):
        return col.value if hasattr(col, "value") else col

    def get_bp_rp(row):
        """Compute bp_rp: use bp_rp if valid, else phot_bp_mean_mag - phot_rp_mean_mag."""
        if "bp_rp" in row.colnames:
            try:
                val = get_value(row["bp_rp"])
                val_float = float(val)
                if np.isfinite(val_float):
                    return val_float
            except (ValueError, TypeError):
                pass
        if "phot_bp_mean_mag" in row.colnames and "phot_rp_mean_mag" in row.colnames:
            try:
                bp = get_value(row["phot_bp_mean_mag"])
                rp = get_value(row["phot_rp_mean_mag"])
                if np.isfinite(float(bp)) and np.isfinite(float(rp)):
                    return float(bp) - float(rp)
            except (ValueError, TypeError):
                pass
        return None

    return [
        (
            int(row["source_id"]),
            float(get_value(row["ra"])),
            float(get_value(row["dec"])),
            float(get_value(row["parallax"])),
            float(get_value(row["phot_g_mean_mag"])),
            get_bp_rp(row),
        )
        for row in table
    ]


def _download_from_gaia(cache_db: Path, chunk_size: int, star_limit=None, existing_count: int = 0):
    """Download data from Gaia API and store in SQLite cache.

    ⚠️ IMPORTANT: This is the ONLY function that queries Gaia API directly.
    All other operations (processing, filtering, plotting) use the SQLite cache.

    Args:
        cache_db: Path to cache database.
        chunk_size: Number of stars per chunk to request from Gaia.
        star_limit: Maximum number of stars to download (None = no limit)
        existing_count: Number of stars already in cache (for progress tracking)
    """
    # Import here to avoid a hard dependency for callers that don't need tqdm/np
    import sys

    import numpy as np

    conn = init_database(cache_db)
    cursor = conn.cursor()

    if star_limit is not None:
        print(f"Downloading Gaia DR3 data (limit: {star_limit:,} stars)...")
    else:
        print("Downloading Gaia DR3 data...")
        print("(Skipping count query - will download until no more results)")

    chunk_num = 0
    offset = 0
    stars_downloaded = 0
    stars_inserted = existing_count  # Track actually inserted (excluding duplicates)

    # Create progress bar
    if HAS_TQDM:
        if star_limit is not None:
            pbar = tqdm(
                total=star_limit,
                initial=existing_count,
                unit="stars",
                desc="Downloading",
                unit_scale=True,
            )
        else:
            pbar = tqdm(
                initial=existing_count,
                unit="stars",
                desc="Downloading",
                unit_scale=True,
            )
    else:
        pbar = None
        if star_limit is not None:
            percent = (stars_inserted / star_limit) * 100 if star_limit > 0 else 0
            print(
                f"Progress: {stars_inserted:,}/{star_limit:,} stars in cache ({percent:.1f}%)"
            )
        else:
            print(f"Progress: {stars_inserted:,} stars in cache")
        sys.stdout.flush()

    try:
        while True:
            if star_limit is not None and stars_inserted >= star_limit:
                print(
                    f"\nReached target of {star_limit:,} stars (inserted: {stars_inserted:,})"
                )
                break

            if HAS_TQDM and pbar is not None:
                pbar.set_description(f"Downloading chunk {chunk_num + 1}")
            else:
                print(
                    f"Chunk {chunk_num + 1}: Downloading... (Total so far: {stars_inserted:,} unique stars)"
                )
                sys.stdout.flush()

            query = f"""
            SELECT TOP {chunk_size} source_id, ra, dec, parallax, phot_g_mean_mag,
                   bp_rp, phot_bp_mean_mag, phot_rp_mean_mag
            FROM gaiadr3.gaia_source
            WHERE phot_g_mean_mag IS NOT NULL AND phot_rp_mean_mag < 16
            ORDER BY random_index
            OFFSET {offset}
            """

            results = None
            last_error = None
            for attempt in range(GAIA_MAX_RETRIES):
                try:
                    # Use async to avoid 408: sync jobs timeout on ORDER BY random_index + OFFSET.
                    job = Gaia.launch_job_async(query)
                    results = job.get_results()
                    last_error = None
                    break
                except Exception as e:
                    last_error = e
                    is_408 = "408" in str(e) or (getattr(e, "response", None) and getattr(getattr(e, "response", None), "status_code", None) == 408)
                    if is_408 and attempt < GAIA_MAX_RETRIES - 1:
                        wait_s = (attempt + 1) * 5
                        print(f"  Gaia job timeout (408), retrying in {wait_s}s (attempt {attempt + 1}/{GAIA_MAX_RETRIES})...")
                        time.sleep(wait_s)
                    else:
                        print(f"  Error during Gaia query: {e}")
                        break
            if last_error is not None:
                break

            if len(results) == 0:
                print("  No more results from Gaia.")
                break

            data_to_insert = _extract_gaia_data_from_table(results)
            if not data_to_insert:
                print("  No valid rows in this chunk.")
                break

            cursor.executemany(
                """
                INSERT OR IGNORE INTO gaia_source
                (source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                data_to_insert,
            )
            conn.commit()

            chunk_inserted = cursor.rowcount
            duplicates = len(data_to_insert) - chunk_inserted
            stars_inserted += chunk_inserted
            stars_downloaded += len(data_to_insert)

            if HAS_TQDM and pbar is not None:
                if star_limit is not None:
                    remaining = max(0, star_limit - pbar.n)
                    if remaining > 0:
                        pbar.update(min(chunk_inserted, remaining))
                else:
                    pbar.update(chunk_inserted)
                pbar.set_postfix(
                    {
                        "chunk": chunk_num + 1,
                        "new": f"{chunk_inserted:,}",
                        "dups": f"{duplicates:,}",
                        "total": f"{stars_inserted:,}",
                    }
                )
            else:
                if star_limit is not None:
                    percent = (
                        stars_inserted / star_limit * 100 if star_limit > 0 else 0
                    )
                    dup_msg = (
                        f", {duplicates:,} duplicates" if duplicates > 0 else ""
                    )
                    print(
                        f"  ✓ Chunk {chunk_num + 1}: Inserted {chunk_inserted:,} new stars{dup_msg} "
                        f"(Total: {stars_inserted:,}/{star_limit:,} stars, {percent:.1f}%)"
                    )
                else:
                    dup_msg = (
                        f", {duplicates:,} duplicates" if duplicates > 0 else ""
                    )
                    print(
                        f"  ✓ Chunk {chunk_num + 1}: Inserted {chunk_inserted:,} new stars{dup_msg} "
                        f"(Total: {stars_inserted:,} stars)"
                    )
                sys.stdout.flush()

            chunk_num += 1
            offset += chunk_size

    finally:
        if HAS_TQDM and pbar is not None:
            pbar.close()
        conn.close()

    # Return the total number of stars now present in the cache. This matches
    # the contract expected by callers such as ensure_cache_populated.
    return stars_inserted