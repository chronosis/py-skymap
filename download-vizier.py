#!/usr/bin/env python3
"""
Bulk download a VizieR catalog into the local SQLite cache.

For Gaia DR3 (I/355/gaiadr3), rows are merged into gaia_source so skymap-gen
can use them like the Gaia Archive or CSV downloads. Other VizieR catalogues
are not yet merged; only Gaia DR3 is supported for gaia_source.

Usage:
  poetry run python download-vizier.py
  poetry run python download-vizier.py I/355/gaiadr3 100000
  poetry run python download-vizier.py I/355/gaiadr3 --limit 50000
  poetry run python download-vizier.py I/355/gaiadr3 --all   # full catalog (very large)
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.constants import CACHE_DB
from lib.vizier_client import GAIA_VIZIER_CATALOG, download_vizier_catalog


def main():
    parser = argparse.ArgumentParser(
        description="Bulk download a VizieR catalog into the gaia_cache SQLite database."
    )
    parser.add_argument(
        "catalog",
        nargs="?",
        default=GAIA_VIZIER_CATALOG,
        help=f"VizieR catalog ID (default: {GAIA_VIZIER_CATALOG})",
    )
    parser.add_argument(
        "limit",
        nargs="?",
        type=int,
        default=None,
        help="Max rows to download (default: 50000). Ignored if --all is set.",
    )
    parser.add_argument(
        "--limit",
        "-n",
        type=int,
        default=None,
        dest="limit_opt",
        help="Max rows to download (overrides positional limit)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Download full catalog (Vizier.ROW_LIMIT = -1). Can be very large for Gaia DR3.",
    )
    parser.add_argument(
        "--no-merge",
        action="store_true",
        help="Do not merge into gaia_source (only fetch; no DB write for now)",
    )
    args = parser.parse_args()

    catalog_id = args.catalog.strip()
    if args.all:
        row_limit = -1
        if "355" in catalog_id and "gaiadr3" in catalog_id.lower():
            print("Warning: Gaia DR3 has ~1.8e9 rows; --all can take a long time.", file=sys.stderr)
    else:
        row_limit = args.limit_opt if args.limit_opt is not None else args.limit
        if row_limit is None:
            row_limit = 50_000

    print(f"Downloading VizieR catalog {catalog_id} (limit={row_limit}) into {CACHE_DB}...")
    try:
        n = download_vizier_catalog(
            catalog_id,
            CACHE_DB,
            row_limit=row_limit,
            merge_into_gaia=not args.no_merge,
        )
        print(f"Inserted {n:,} rows into gaia_source.")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
