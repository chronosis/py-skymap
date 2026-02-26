#!/usr/bin/env python3
"""
One-off script to drop the star_positions_3d table from the cache database.

This table stored target-relative coordinates and is no longer used. Run this
to clean up existing databases after migrating to Earth-centric star_cartesian_earth.

Usage:
    python drop_star_positions_3d.py [--db PATH]

If --db is omitted, uses CACHE_DB from lib.constants (gaia_cache/gaia_cache.db).
"""
import argparse
import sqlite3
import sys

from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Drop star_positions_3d table from cache DB")
    parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help="Path to cache database (default: gaia_cache/gaia_cache.db)",
    )
    args = parser.parse_args()

    if args.db is not None:
        db_path = args.db.resolve()
    else:
        from lib.constants import CACHE_DB
        db_path = Path.cwd() / CACHE_DB

    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        sys.exit(1)

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Check if table exists
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='star_positions_3d'"
    )
    if cursor.fetchone() is None:
        print(f"Table star_positions_3d does not exist in {db_path}")
        conn.close()
        return 0

    # Get row count before drop (informational)
    cursor.execute("SELECT COUNT(*) FROM star_positions_3d")
    count = cursor.fetchone()[0]
    print(f"Dropping star_positions_3d ({count:,} rows) from {db_path}...")

    cursor.execute("DROP TABLE star_positions_3d")
    conn.commit()
    conn.close()

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
