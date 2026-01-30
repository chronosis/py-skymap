#!/usr/bin/env python3
"""
Query GAIA data endpoint and output each row as a JSON blob to the console.

Configure the query in the QUERY variable below.
"""

import json
import sys
from astroquery.gaia import Gaia
from astropy.table import Table

from lib.constants import GAIA_MAX_RETRIES

# Configure your query here
QUERY = """
SELECT TOP 1 *
FROM gaiadr3.gaia_source
ORDER BY phot_g_mean_mag DESC
"""

def table_row_to_dict(row):
    """Convert an astropy Table row to a dictionary, handling units properly."""
    row_dict = {}
    for colname in row.colnames:
        value = row[colname]
        # Handle astropy Quantity objects (values with units)
        if hasattr(value, 'value'):
            row_dict[colname] = value.value
        elif hasattr(value, 'item'):
            # Handle numpy scalars
            row_dict[colname] = value.item()
        else:
            row_dict[colname] = value
    return row_dict

def query_gaia_and_output_json(query, max_retries=GAIA_MAX_RETRIES):
    """
    Execute a GAIA query and output each row as a JSON blob to stdout.
    
    Parameters:
    -----------
    query : str
        The ADQL query to execute
    max_retries : int
        Maximum number of retry attempts for timeouts
    """
    print(f"Executing GAIA query...", file=sys.stderr)
    print(f"Query:\n{query}\n", file=sys.stderr)
    
    retry_count = 0
    job = None
    results = None
    
    while retry_count < max_retries:
        try:
            job = Gaia.launch_job(query)
            results = job.get_results()
            break  # Success, exit retry loop
        except Exception as e:
            retry_count += 1
            if "timeout" in str(e).lower() or "408" in str(e) or retry_count >= max_retries:
                if retry_count >= max_retries:
                    print(f"Error: Failed after {max_retries} retries. {e}", file=sys.stderr)
                    sys.exit(1)
                else:
                    wait_time = retry_count * 5
                    print(f"Timeout occurred, retrying in {wait_time}s... (attempt {retry_count}/{max_retries})", file=sys.stderr)
                    import time
                    time.sleep(wait_time)
            else:
                raise  # Re-raise if it's not a timeout
    
    if results is None or len(results) == 0:
        print("No results returned from query.", file=sys.stderr)
        return
    
    print(f"Found {len(results)} rows. Outputting JSON...\n", file=sys.stderr)
    
    # Output each row as a JSON blob
    for row in results:
        row_dict = table_row_to_dict(row)
        # Output as a single-line JSON blob
        print(json.dumps(row_dict))
    
    print(f"\nCompleted. Output {len(results)} rows.", file=sys.stderr)

if __name__ == "__main__":
    query_gaia_and_output_json(QUERY)
