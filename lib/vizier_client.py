"""
Bulk download from VizieR and optionally merge into the gaia_source cache.

VizieR hosts many published catalogues (including Gaia DR3 as I/355/gaiadr3).
Use this to bulk-download star data when you prefer VizieR over the Gaia Archive
or want to pull other VizieR catalogues into the same SQLite cache.
"""

from pathlib import Path

from astropy.table import Table

from .sqlite_helper import init_database

# Gaia DR3 on VizieR: map VizieR column names to gaia_source schema.
# VizieR may use Source, RA_ICRS, DE_ICRS, Plx, Gmag, BP-RP (or BP-RP_).
GAIA_VIZIER_CATALOG = "I/355/gaiadr3"
GAIA_VIZIER_COLUMN_ALIASES = {
    "source_id": ["Source", "source_id", "SOURCE"],
    "ra": ["RA_ICRS", "ra", "RA"],
    "dec": ["DE_ICRS", "dec", "DEC", "DE_ICRS"],
    "parallax": ["Plx", "parallax", "PLX"],
    "phot_g_mean_mag": ["Gmag", "phot_g_mean_mag", "Gmag_"],
    "bp_rp": ["BP-RP", "BP-RP_", "bp_rp", "BP_RP"],
}


def _get_col(table: Table, aliases: list[str]):
    """Return first existing column name from table for the given aliases."""
    for name in aliases:
        if name in table.colnames:
            return name
    return None


def _extract_gaia_like_rows(table: Table) -> list[tuple]:
    """Convert a VizieR table (Gaia DR3 or similar) to rows for gaia_source.

    Returns list of (source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp).
    Rows with missing required fields are skipped.
    """
    col_src = _get_col(table, GAIA_VIZIER_COLUMN_ALIASES["source_id"])
    col_ra = _get_col(table, GAIA_VIZIER_COLUMN_ALIASES["ra"])
    col_dec = _get_col(table, GAIA_VIZIER_COLUMN_ALIASES["dec"])
    col_plx = _get_col(table, GAIA_VIZIER_COLUMN_ALIASES["parallax"])
    col_gmag = _get_col(table, GAIA_VIZIER_COLUMN_ALIASES["phot_g_mean_mag"])
    col_bprp = _get_col(table, GAIA_VIZIER_COLUMN_ALIASES["bp_rp"])

    if not all([col_src, col_ra, col_dec, col_plx, col_gmag]):
        return []

    def _val(row, col, default=None):
        if col is None:
            return default
        v = row[col]
        if hasattr(v, "value"):
            v = v.value
        try:
            if hasattr(v, "filled"):
                v = v.filled(float("nan"))
            if v is None or (isinstance(v, float) and (v != v or v == float("inf"))):
                return default
            return float(v)
        except (TypeError, ValueError):
            return default

    rows = []
    for row in table:
        try:
            sid = int(_val(row, col_src, 0))
            ra = _val(row, col_ra)
            dec = _val(row, col_dec)
            plx = _val(row, col_plx)
            gmag = _val(row, col_gmag)
            bprp = _val(row, col_bprp)
        except (TypeError, ValueError):
            continue
        if ra is None or dec is None or plx is None or gmag is None:
            continue
        if not (sid and abs(ra) <= 360 and abs(dec) <= 90):
            continue
        rows.append((sid, float(ra), float(dec), float(plx), float(gmag), bprp))
    return rows


def download_vizier_catalog(
    catalog_id: str,
    cache_db: Path,
    *,
    row_limit: int | None = 50_000,
    merge_into_gaia: bool = True,
) -> int:
    """
    Download a VizieR catalog and, if it is Gaia DR3 (I/355/gaiadr3), merge into
    gaia_source in the SQLite cache.

    Args:
        catalog_id: VizieR catalog identifier (e.g. "I/355/gaiadr3").
        cache_db: Path to the cache database.
        row_limit: Max rows to download (None = use Vizier default; -1 = unlimited).
        merge_into_gaia: If True and catalog is Gaia DR3, insert into gaia_source.

    Returns:
        Number of rows inserted into gaia_source (0 if not merged or no rows).
    """
    from astroquery.vizier import Vizier

    Vizier.ROW_LIMIT = row_limit if row_limit is not None else 50
    if row_limit == -1:
        Vizier.ROW_LIMIT = -1

    tables = Vizier.get_catalogs(catalog_id)
    if not tables:
        return 0

    # TableList / dict: use first table
    if hasattr(tables, "keys"):
        first_key = next(iter(tables.keys()), None)
        table = tables[first_key] if first_key is not None else None
    elif hasattr(tables, "__getitem__"):
        table = tables[0]
    else:
        table = tables
    if not isinstance(table, Table) or len(table) == 0:
        return 0

    catalog_norm = catalog_id.strip().upper().replace(" ", "").replace("\\", "/")
    is_gaia_dr3 = "I/355/GAIADR3" in catalog_norm or catalog_norm == "I/355/GAIADR3"
    if merge_into_gaia and is_gaia_dr3:
        rows = _extract_gaia_like_rows(table)
        if not rows:
            return 0
        conn = init_database(cache_db)
        cursor = conn.cursor()
        cursor.executemany(
            """
            INSERT OR IGNORE INTO gaia_source
            (source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            rows,
        )
        n = cursor.rowcount
        conn.commit()
        conn.close()
        return n
    return 0
