#!/usr/bin/env python3
"""
One-time migration: copy all data from the SQLite cache into PostgreSQL.

Usage:
    python migrate_sqlite_to_pg.py [--sqlite PATH] [--pg-dsn DSN] [--chunk-size N]

Defaults:
    --sqlite   gaia_cache/gaia_cache.db  (CACHE_DB from lib.constants)
    --pg-dsn   $PG_DSN env var           (PG_DSN from lib.constants)
    --chunk-size 50000

The script is safely re-runnable: it uses ON CONFLICT DO NOTHING so rows
already present in PostgreSQL are skipped.
"""

import argparse
import sqlite3
import sys
import time

import psycopg2
import psycopg2.extras

from lib.constants import CACHE_DB, PG_DSN
from lib.pg_helper import init_database as pg_init_database


TABLES = [
    {
        "name": "gaia_source",
        "columns": "source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp, source",
        "conflict_col": "source_id",
    },
    {
        "name": "star_cartesian_earth",
        "columns": "source_id, x_pc, y_pc, z_pc, source",
        "conflict_col": "source_id",
    },
    {
        "name": "simbad_cache",
        "columns": "main_id, ra, dec, parallax_mas, vmag, otype, cached_at",
        "conflict_col": "main_id",
    },
    {
        "name": "simbad_cache_aliases",
        "columns": "alias, main_id",
        "conflict_col": "alias",
    },
]


def _table_exists_sqlite(sqlite_cur, table_name: str) -> bool:
    sqlite_cur.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (table_name,),
    )
    return sqlite_cur.fetchone() is not None


def _count_sqlite(sqlite_cur, table_name: str) -> int:
    sqlite_cur.execute(f"SELECT COUNT(*) FROM {table_name}")
    return sqlite_cur.fetchone()[0]


def migrate_table(
    sqlite_conn,
    pg_conn,
    table_name: str,
    columns: str,
    conflict_col: str,
    chunk_size: int,
) -> int:
    sqlite_cur = sqlite_conn.cursor()

    if not _table_exists_sqlite(sqlite_cur, table_name):
        print(f"  {table_name}: table does not exist in SQLite, skipping.")
        return 0

    total = _count_sqlite(sqlite_cur, table_name)
    if total == 0:
        print(f"  {table_name}: 0 rows, skipping.")
        return 0

    print(f"  {table_name}: {total:,} rows to migrate...")

    col_list = [c.strip() for c in columns.split(",")]
    n_cols = len(col_list)
    placeholders = ", ".join(["%s"] * n_cols)

    insert_sql = (
        f"INSERT INTO {table_name} ({columns}) VALUES %s "
        f"ON CONFLICT ({conflict_col}) DO NOTHING"
    )

    pg_cur = pg_conn.cursor()
    offset = 0
    total_inserted = 0
    t0 = time.time()

    while offset < total:
        sqlite_cur.execute(
            f"SELECT {columns} FROM {table_name} LIMIT {chunk_size} OFFSET {offset}"
        )
        rows = sqlite_cur.fetchall()
        if not rows:
            break

        psycopg2.extras.execute_values(
            pg_cur, insert_sql, rows, page_size=5000
        )
        pg_conn.commit()
        inserted = pg_cur.rowcount
        total_inserted += inserted
        offset += len(rows)

        elapsed = time.time() - t0
        pct = offset / total * 100
        rate = offset / elapsed if elapsed > 0 else 0
        print(
            f"    {offset:,}/{total:,} read ({pct:.1f}%), "
            f"{total_inserted:,} new, "
            f"{rate:,.0f} rows/sec"
        )

    elapsed = time.time() - t0
    print(
        f"  {table_name}: done. {total_inserted:,} new rows inserted in {elapsed:.1f}s"
    )
    return total_inserted


def main():
    parser = argparse.ArgumentParser(description="Migrate SQLite cache to PostgreSQL")
    parser.add_argument(
        "--sqlite",
        type=str,
        default=str(CACHE_DB),
        help=f"Path to SQLite database (default: {CACHE_DB})",
    )
    parser.add_argument(
        "--pg-dsn",
        type=str,
        default=PG_DSN,
        help="PostgreSQL DSN (default: PG_DSN from .env)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=50_000,
        help="Rows per batch (default: 50000)",
    )
    args = parser.parse_args()

    from pathlib import Path

    sqlite_path = Path(args.sqlite)
    if not sqlite_path.exists():
        print(f"SQLite database not found: {sqlite_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Source:  {sqlite_path} ({sqlite_path.stat().st_size / 1e9:.2f} GB)")
    print(f"Target:  {args.pg_dsn.split('@')[-1] if '@' in args.pg_dsn else args.pg_dsn}")
    print(f"Chunk:   {args.chunk_size:,} rows\n")

    print("Ensuring PostgreSQL schema exists...")
    pg_conn = pg_init_database(args.pg_dsn)
    pg_conn.close()

    sqlite_conn = sqlite3.connect(str(sqlite_path))
    pg_conn = psycopg2.connect(args.pg_dsn)

    t_start = time.time()
    grand_total = 0

    for tbl in TABLES:
        n = migrate_table(
            sqlite_conn,
            pg_conn,
            tbl["name"],
            tbl["columns"],
            tbl["conflict_col"],
            args.chunk_size,
        )
        grand_total += n

    sqlite_conn.close()
    pg_conn.close()

    elapsed = time.time() - t_start
    print(f"\nMigration complete: {grand_total:,} new rows in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
