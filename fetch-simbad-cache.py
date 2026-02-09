#!/usr/bin/env python3
"""
Populate the Simbad cache from a list of object names (default: bright stars from star_data).

Run once to fill the cache; skymap-gen and other code can then use Simbad-sourced
positions/magnitudes/parallax to flesh out gaps. Simbad is rate-limited; this script
batches requests and adds a short delay between batches.

Usage:
  poetry run python fetch-simbad-cache.py
  poetry run python fetch-simbad-cache.py "M1,M2,NGC 1976"
  poetry run python fetch-simbad-cache.py --names "Sirius,Aldebaran,Vega"
"""

import argparse
import sys
from pathlib import Path

# Add project root so lib is importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.constants import CACHE_DB
from lib.simbad_client import query_simbad_and_cache
from lib.star_data import get_bright_stars


def main():
    parser = argparse.ArgumentParser(
        description="Fetch Simbad data for named objects and cache in SQLite."
    )
    parser.add_argument(
        "names",
        nargs="?",
        default=None,
        help="Comma-separated object names (default: use built-in bright star list)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50,
        help="Number of objects per Simbad batch query (default: 50)",
    )
    args = parser.parse_args()

    if args.names:
        identifiers = [s.strip() for s in args.names.split(",") if s.strip()]
    else:
        identifiers = [s["name"] for s in get_bright_stars()]

    if not identifiers:
        print("No identifiers to fetch.")
        return 0

    print(f"Fetching Simbad data for {len(identifiers)} object(s) into {CACHE_DB}...")
    total = 0
    for i in range(0, len(identifiers), args.batch_size):
        batch = identifiers[i : i + args.batch_size]
        try:
            rows = query_simbad_and_cache(CACHE_DB, batch)
            total += len(rows)
            print(f"  Batch {i // args.batch_size + 1}: cached {len(rows)}/{len(batch)}")
        except Exception as e:
            print(f"  Batch {i // args.batch_size + 1} failed: {e}", file=sys.stderr)
    print(f"Done. Total cached: {total} object(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
