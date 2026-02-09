"""
Query Simbad and cache results in SQLite, similar to Gaia cache.

Use this to flesh out gaps: e.g. objects not in Gaia, or to refresh positions/magnitudes
for known objects. Simbad is rate-limited (recommended 5–10 queries/sec); batch
queries via query_objects() count as one query.
"""

import time
from pathlib import Path

from astropy.table import Table

from .sqlite_helper import init_database, get_simbad_from_cache, put_simbad_in_cache

# Rate limit: seconds between Simbad requests (batch counts as one).
SIMBAD_QUERY_DELAY_SEC = 0.2


def _row_to_dict(table: Table, row) -> dict | None:
    """Convert one row of a Simbad result table to our cache dict.

    Handles both uppercase and lowercase column names.
    """
    def _get_val(val):
        if val is None:
            return None
        if hasattr(val, "value"):
            val = val.value
        try:
            if hasattr(val, "filled"):
                val = val.filled()
            if val is None or (isinstance(val, float) and (val != val or val == float("inf"))):
                return None
            return val
        except Exception:
            return None

    def col(name: str):
        for key in table.colnames:
            if key.upper() == name.upper():
                val = _get_val(row[key])
                if val is None:
                    return None
                try:
                    return float(val)
                except (TypeError, ValueError):
                    return None
        return None

    def col_str(name: str):
        for key in table.colnames:
            if key.upper() == name.upper():
                val = _get_val(row[key])
                return str(val).strip() if val is not None else None
        return None

    main_id = col_str("MAIN_ID") or col_str("main_id")
    if not main_id:
        return None
    ra = col("RA") if col("RA") is not None else col("ra")
    dec = col("DEC") if col("DEC") is not None else col("dec")
    if ra is None or dec is None:
        return None

    out = {
        "main_id": main_id,
        "ra": float(ra),
        "dec": float(dec),
        "parallax_mas": col("PLX_VALUE") or col("PLX") or col("PARALLAX"),
        "vmag": col("V") or col("VMAG"),
        "otype": col_str("OTYPE"),
    }
    return out


def query_simbad_and_cache(
    cache_db: Path,
    identifiers: list[str],
    *,
    delay_sec: float = SIMBAD_QUERY_DELAY_SEC,
) -> list[dict]:
    """
    Query Simbad for the given object names, insert results into the Simbad cache,
    and return a list of dicts (main_id, ra, dec, parallax_mas, vmag, otype) for each
    successfully resolved object.

    Uses astroquery.simbad.Simbad with parallax and magnitude fields. Batch request
    is used (one query per call) to respect rate limits. Results are cached in
    simbad_cache table.
    """
    from astroquery.simbad import Simbad

    if not identifiers:
        return []

    # Ensure cache table exists
    init_database(cache_db)

    # Add extra fields so we get parallax and magnitude (Simbad uses 'plx' and 'V')
    Simbad.add_votable_fields("plx", "V", "otype")

    # Batch query (one request for all names)
    identifiers_clean = [s.strip() for s in identifiers if s and s.strip()]
    if not identifiers_clean:
        return []

    try:
        result = Simbad.query_objects(identifiers_clean)
    except Exception as e:
        raise RuntimeError(f"Simbad query failed: {e}") from e

    if delay_sec > 0:
        time.sleep(delay_sec)

    if result is None or len(result) == 0:
        return []

    rows = []
    for row in result:
        d = _row_to_dict(result, row)
        if d:
            rows.append(d)

    if rows:
        # Only store aliases when result order matches request order (same count).
        requested = identifiers_clean if len(rows) == len(identifiers_clean) else None
        put_simbad_in_cache(cache_db, rows, requested_identifiers=requested)
    return rows


def get_from_simbad_or_cache(cache_db: Path, identifier: str) -> dict | None:
    """
    Return Simbad data for one object: from cache if present, otherwise query
    Simbad, cache the result, and return it. Returns None if not found or on error.
    """
    init_database(cache_db)
    cached = get_simbad_from_cache(cache_db, identifier)
    if cached is not None:
        return cached

    from astroquery.simbad import Simbad

    Simbad.add_votable_fields("plx", "V", "otype")
    try:
        result = Simbad.query_object(identifier.strip())
    except Exception:
        return None

    if result is None or len(result) == 0:
        return None

    time.sleep(SIMBAD_QUERY_DELAY_SEC)

    d = _row_to_dict(result, result[0])
    if d:
        put_simbad_in_cache(
            cache_db,
            [d],
            requested_identifiers=[identifier.strip()],
        )
        return d
    return None


def simbad_entry_to_bright_star_format(entry: dict) -> dict:
    """
    Convert a simbad_cache entry (or get_from_simbad_or_cache result) to the
    same shape as get_bright_stars() items: name, ra_deg, dec_deg, apparent_mag, parallax_mas.
    """
    return {
        "name": entry["main_id"],
        "ra_deg": entry["ra"],
        "dec_deg": entry["dec"],
        "apparent_mag": entry.get("vmag"),
        "parallax_mas": entry.get("parallax_mas"),
    }
