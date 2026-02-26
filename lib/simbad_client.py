"""
Query Simbad and cache results in SQLite, similar to Gaia cache.

Use this to flesh out gaps: e.g. objects not in Gaia, or to refresh positions/magnitudes
for known objects. Simbad is rate-limited (recommended 5–10 queries/sec); batch
queries via query_objects() count as one query.
"""

import hashlib
import time
from pathlib import Path

from astropy.table import Table

from .constants import ASSUMED_BACKGROUND_DISTANCE_PC
from .sqlite_helper import (
    init_database,
    get_simbad_from_cache,
    put_simbad_in_cache,
    ra_dec_distance_to_cartesian,
    insert_star_cartesian_earth_batch,
)

# Rate limit: seconds between Simbad requests (batch counts as one).
SIMBAD_QUERY_DELAY_SEC = 0.2

# Explicit TAP base URL for Simbad (used when astroquery's Simbad has no/invalid TAP URL).
SIMBAD_TAP_URL = "https://simbad.cds.unistra.fr/simbad/sim-tap"


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


def _simbad_source_id(main_id: str) -> int:
    """Stable negative source_id for Simbad rows (avoids clashing with Gaia IDs)."""
    h = int(hashlib.md5(main_id.encode("utf-8")).hexdigest()[:8], 16)
    return -((h % (2**30 - 1)) + 1)


def _tap_row_to_dict(table: Table, row) -> dict | None:
    """Convert one row of a Simbad TAP result table to our cache dict.

    TAP may return columns like main_id, ra, dec, plx (or plx_value); reuses
    the same case-insensitive lookup as _row_to_dict.
    """
    return _row_to_dict(table, row)


def query_simbad_bulk_tap(
    cache_db: Path,
    row_limit: int,
    *,
    async_job: bool = False,
    timeout_sec: int = 300,
    otype_filter: str | None = "Star",
) -> list[dict]:
    """
    Query Simbad via TAP for an arbitrary number of stars in one bulk request.

    No pre-defined star list is required; Simbad returns up to row_limit rows
    from its catalog (objects with valid ra/dec; optionally filtered by type).
    Result rows are returned in the same dict shape as query_simbad_and_cache
    (main_id, ra, dec, parallax_mas, vmag, otype).

    Args:
        cache_db: Path to cache database (used to ensure schema exists).
        row_limit: Maximum number of rows to request (TOP N in ADQL).
        async_job: If True, run the TAP job asynchronously. Defaults to False to
            avoid issues with some TAP async endpoints.
        timeout_sec: Timeout in seconds when async_job is True.
        otype_filter: Restrict to objects whose otype starts with this string
            (e.g. "Star"); use None to skip type filter.

    Returns:
        List of dicts with keys main_id, ra, dec, parallax_mas, vmag, otype
        (vmag/parallax_mas may be None if not in TAP response).
    """
    from astroquery.simbad import Simbad

    init_database(cache_db)

    # Simbad TAP "basic" table: main_id, ra, dec, plx (parallax in mas). otype for object type.
    # Filter to valid coordinates; optional otype to focus on stars.
    where_parts = ["ra IS NOT NULL", "dec IS NOT NULL"]
    if otype_filter:
        where_parts.append(f"otype LIKE '{otype_filter}%'")
    where_clause = " AND ".join(where_parts)

    adql = (
        f"SELECT TOP {int(row_limit)} main_id, ra, dec, plx, otype "
        f"FROM basic WHERE {where_clause}"
    )

    time.sleep(SIMBAD_QUERY_DELAY_SEC)

    table = None
    try:
        result = Simbad.query_tap(
            adql,
            async_job=async_job,
            timeout=timeout_sec if async_job else None,
        )
        if result is not None:
            table = result.to_table() if hasattr(result, "to_table") else result
    except Exception as e:
        err_msg = str(e)
        # Astroquery may leave TAP URL unset (None), leading to "Bad URI syntax: 'None'"
        if "None" in err_msg or "URI" in err_msg or "uri" in err_msg:
            try:
                from pyvo.dal import TAPService
                tap = TAPService(SIMBAD_TAP_URL)
                result = tap.search(adql, language="ADQL", maxrec=row_limit)
                if result is not None and hasattr(result, "to_table"):
                    table = result.to_table()
            except Exception as pyvo_err:
                raise RuntimeError(
                    f"Simbad TAP query failed (astroquery: {e}; pyvo fallback: {pyvo_err})"
                ) from pyvo_err
        else:
            raise RuntimeError(f"Simbad TAP query failed: {e}") from e

    if table is None or not isinstance(table, Table) or len(table) == 0:
        return []

    # _row_to_dict uses case-insensitive col("PLX"), col("PLX_VALUE"), col("PARALLAX")
    rows = []
    for row in table:
        d = _tap_row_to_dict(table, row)
        if d:
            rows.append(d)
    return rows


def ingest_simbad_bulk_into_gaia_source(
    cache_db: Path,
    row_limit: int,
    *,
    async_job: bool = False,
    timeout_sec: int = 300,
    cache_simbad: bool = True,
) -> int:
    """
    Fetch an arbitrary number of stars from Simbad via TAP and insert them into
    gaia_source with source='simbad'. Does not require a pre-defined star list.

    Args:
        cache_db: Path to cache database.
        row_limit: Max number of stars to request from Simbad TAP.
        async_job: Passed to query_simbad_bulk_tap (defaults to False).
        timeout_sec: Passed to query_simbad_bulk_tap.
        cache_simbad: If True, also write results into simbad_cache for lookups.

    Returns:
        Number of new rows inserted into gaia_source.
    """
    rows_simbad = query_simbad_bulk_tap(
        cache_db,
        row_limit,
        async_job=async_job,
        timeout_sec=timeout_sec,
    )
    if not rows_simbad:
        return 0

    if cache_simbad:
        put_simbad_in_cache(cache_db, rows_simbad, requested_identifiers=None)

    gaia_rows = []
    for d in rows_simbad:
        main_id = d.get("main_id")
        if not main_id:
            continue
        ra = d.get("ra")
        dec = d.get("dec")
        if ra is None or dec is None:
            continue
        source_id = _simbad_source_id(main_id)
        parallax = d.get("parallax_mas")
        if parallax is None or (
            isinstance(parallax, float) and (parallax != parallax or parallax <= 0)
        ):
            parallax = 0.0
        mag = d.get("vmag")
        if mag is None or (isinstance(mag, float) and (mag != mag)):
            mag = 99.0
        gaia_rows.append(
            (source_id, float(ra), float(dec), float(parallax), float(mag), None)
        )

    if not gaia_rows:
        return 0

    conn = init_database(cache_db)
    cursor = conn.cursor()
    cursor.executemany(
        """
        INSERT OR IGNORE INTO gaia_source
        (source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp, source)
        VALUES (?, ?, ?, ?, ?, ?, 'simbad')
        """,
        gaia_rows,
    )
    n = cursor.rowcount
    conn.commit()
    conn.close()
    # Populate Earth-centric Cartesian cache
    cartesian_rows = []
    for (sid, ra, dec, plx, gmag, bp_rp) in gaia_rows:
        if plx and plx > 0 and plx == plx:
            d_pc = 1000.0 / plx
        else:
            d_pc = ASSUMED_BACKGROUND_DISTANCE_PC
        x, y, z = ra_dec_distance_to_cartesian(ra, dec, d_pc)
        cartesian_rows.append((sid, x, y, z, "simbad"))
    if cartesian_rows:
        insert_star_cartesian_earth_batch(cache_db, cartesian_rows)
    return n


def ingest_simbad_into_gaia_source(
    cache_db: Path,
    identifiers: list[str],
) -> int:
    """
    Ensure Simbad cache has data for the given identifiers, then merge those
    entries into gaia_source with source='simbad' so they appear in star counts
    and on the map. Uses stable negative source_ids so Simbad rows do not clash
    with Gaia. Missing parallax/magnitude use 0.0 and 99.0 respectively.

    Returns:
        Number of rows inserted into gaia_source (new Simbad-sourced rows).
    """
    if not identifiers:
        return 0
    ids_clean = [s.strip() for s in identifiers if s and s.strip()]
    if not ids_clean:
        return 0

    # Populate Simbad cache for these identifiers (no-op if already cached)
    rows_simbad = query_simbad_and_cache(cache_db, ids_clean)
    if not rows_simbad:
        return 0

    gaia_rows = []
    for d in rows_simbad:
        main_id = d.get("main_id")
        if not main_id:
            continue
        ra = d.get("ra")
        dec = d.get("dec")
        if ra is None or dec is None:
            continue
        source_id = _simbad_source_id(main_id)
        parallax = d.get("parallax_mas")
        if parallax is None or (isinstance(parallax, float) and (parallax != parallax or parallax <= 0)):
            parallax = 0.0
        mag = d.get("vmag")
        if mag is None or (isinstance(mag, float) and (mag != mag)):
            mag = 99.0
        gaia_rows.append((source_id, float(ra), float(dec), float(parallax), float(mag), None))

    if not gaia_rows:
        return 0

    conn = init_database(cache_db)
    cursor = conn.cursor()
    cursor.executemany(
        """
        INSERT OR IGNORE INTO gaia_source
        (source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp, source)
        VALUES (?, ?, ?, ?, ?, ?, 'simbad')
        """,
        gaia_rows,
    )
    n = cursor.rowcount
    conn.commit()
    conn.close()
    # Populate Earth-centric Cartesian cache
    cartesian_rows = []
    for (sid, ra, dec, plx, gmag, bp_rp) in gaia_rows:
        if plx and plx > 0 and plx == plx:
            d_pc = 1000.0 / plx
        else:
            d_pc = ASSUMED_BACKGROUND_DISTANCE_PC
        x, y, z = ra_dec_distance_to_cartesian(ra, dec, d_pc)
        cartesian_rows.append((sid, x, y, z, "simbad"))
    if cartesian_rows:
        insert_star_cartesian_earth_batch(cache_db, cartesian_rows)
    return n
