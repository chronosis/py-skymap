#!/usr/bin/env python3
"""
Download Gaia DR3 CSV.gz files from the ESA CDN.
Limits download to approximately the specified number of stars.
"""

import os
import sys
import time
import argparse
import requests
from pathlib import Path
from urllib.parse import urljoin, urlparse
from bs4 import BeautifulSoup

from lib.constants import (
    GAIA_MAX_RETRIES,
    GAIA_HTTP_TIMEOUT_SECONDS,
    GAIA_HTTP_HEAD_TIMEOUT_SECONDS,
)

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(iterable, **kwargs):
        return iterable

BASE_URL = "https://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/"
DOWNLOAD_DIR = Path("gaia_csv")
CHUNK_SIZE = 8192  # Download chunk size in bytes
DEFAULT_STAR_LIMIT = 50000  # Default number of stars to download

def estimate_stars_from_filename(filename):
    """Estimate number of stars from filename range (e.g., GaiaSource_000000-003111.csv.gz)."""
    try:
        # Extract range from filename like "GaiaSource_000000-003111.csv.gz"
        base = filename.replace('.csv.gz', '').replace('GaiaSource_', '')
        if '-' in base:
            start, end = base.split('-')
            # Estimate: not all IDs in range exist, use ~80% as conservative estimate
            return int((int(end) - int(start)) * 0.8)
    except:
        pass
    # Fallback: estimate based on typical file size (~200MB ≈ 200k stars)
    return 200000

def get_file_list():
    """Parse the HTML directory listing to get all CSV.gz file URLs."""
    print(f"Fetching file list from {BASE_URL}...")
    response = requests.get(BASE_URL)
    response.raise_for_status()
    
    soup = BeautifulSoup(response.text, 'html.parser')
    files = []
    
    # Find all links to CSV.gz files
    for link in soup.find_all('a'):
        href = link.get('href', '')
        if href.endswith('.csv.gz'):
            # Decode URL-encoded filenames (e.g., %5F becomes _)
            filename = href.replace('%5F', '_')
            file_url = urljoin(BASE_URL, href)
            files.append((filename, file_url))
    
    files.sort()  # Sort by filename
    print(f"Found {len(files)} CSV.gz files")
    return files

def download_file(url, filepath, max_retries=GAIA_MAX_RETRIES):
    """Download a single file with retry logic."""
    for attempt in range(max_retries):
        try:
            response = requests.get(url, stream=True, timeout=GAIA_HTTP_TIMEOUT_SECONDS)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            # Check if file already exists and is complete
            if filepath.exists():
                existing_size = filepath.stat().st_size
                if existing_size == total_size and total_size > 0:
                    return True  # File already downloaded
            
            # Download with progress
            downloaded = 0
            with open(filepath, 'wb') as f:
                for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
            
            return True
            
        except Exception as e:
            if attempt < max_retries - 1:
                wait_time = (attempt + 1) * 5
                print(f"  ⚠ Error downloading {filepath.name}, retrying in {wait_time}s... (attempt {attempt + 1}/{max_retries})")
                sys.stdout.flush()
                time.sleep(wait_time)
            else:
                print(f"  ✗ Failed to download {filepath.name} after {max_retries} attempts: {e}")
                sys.stdout.flush()
                return False
    
    return False

def main():
    """Main download function."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Download Gaia DR3 CSV.gz files from ESA CDN',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  {sys.argv[0]}              # Download ~{DEFAULT_STAR_LIMIT:,} stars (default)
  {sys.argv[0]} 100000       # Download ~100,000 stars
  {sys.argv[0]} 500000       # Download ~500,000 stars
        """
    )
    parser.add_argument(
        'stars',
        type=int,
        nargs='?',
        default=DEFAULT_STAR_LIMIT,
        help=f'Number of stars to download (default: {DEFAULT_STAR_LIMIT:,})'
    )
    args = parser.parse_args()
    
    star_limit = args.stars
    print(f"Target: Download approximately {star_limit:,} stars")
    
    DOWNLOAD_DIR.mkdir(exist_ok=True)
    
    # Get list of files to download
    files = get_file_list()
    
    if not files:
        print("No CSV.gz files found!")
        return
    
    # Estimate how many files we'll need
    estimated_stars_per_file = estimate_stars_from_filename(files[0][0])
    estimated_files_needed = max(1, int(star_limit / estimated_stars_per_file) + 1)
    files_to_process = files[:min(estimated_files_needed + 5, len(files))]  # Add buffer
    
    print(f"Estimated files needed: ~{estimated_files_needed} (processing {len(files_to_process)} files)")
    print(f"\nDownloading files to {DOWNLOAD_DIR}...")
    
    # Create progress bar
    if HAS_TQDM:
        pbar = tqdm(total=star_limit, desc='Downloading', unit='stars', unit_scale=True)
    else:
        pbar = None
        print(f"Progress: 0/{star_limit:,} stars (0.0%)")
        sys.stdout.flush()
    
    downloaded = 0
    failed = 0
    skipped = 0
    total_stars_estimated = 0
    
    try:
        for i, (filename, url) in enumerate(files_to_process):
            # Check if we've reached the star limit
            if total_stars_estimated >= star_limit:
                print(f"\nReached target of {star_limit:,} stars (estimated: {total_stars_estimated:,})")
                break
            
            filepath = DOWNLOAD_DIR / filename
            estimated_stars = estimate_stars_from_filename(filename)
            
            if HAS_TQDM:
                pbar.set_description(f"Downloading {filename}")
            else:
                percent = (total_stars_estimated / star_limit) * 100 if star_limit > 0 else 0
                print(f"File {i+1}/{len(files_to_process)}: {filename}... "
                      f"({total_stars_estimated:,}/{star_limit:,} stars, {percent:.1f}%)")
                sys.stdout.flush()
            
            # Check if file already exists
            file_was_skipped = False
            if filepath.exists():
                # Try to verify it's complete by checking size
                try:
                    response = requests.head(url, timeout=GAIA_HTTP_HEAD_TIMEOUT_SECONDS)
                    expected_size = int(response.headers.get('content-length', 0))
                    actual_size = filepath.stat().st_size
                    if expected_size > 0 and actual_size == expected_size:
                        skipped += 1
                        file_was_skipped = True
                        total_stars_estimated += estimated_stars
                        if HAS_TQDM:
                            # Update progress, but cap at total
                            remaining = max(0, star_limit - pbar.n)
                            if remaining > 0:
                                pbar.update(min(estimated_stars, remaining))
                            pbar.set_postfix({
                                'skipped': skipped, 
                                'downloaded': downloaded, 
                                'failed': failed,
                                'stars': f'{total_stars_estimated:,}'
                            })
                        else:
                            print(f"  ⊘ Skipped (already exists): {filename} (~{estimated_stars:,} stars)")
                            sys.stdout.flush()
                        continue
                except:
                    pass  # If HEAD request fails, try downloading anyway
            
            # Download the file
            success = download_file(url, filepath)
            
            if success:
                downloaded += 1
                total_stars_estimated += estimated_stars
                if HAS_TQDM:
                    # Update progress, but cap at total
                    remaining = max(0, star_limit - pbar.n)
                    if remaining > 0:
                        pbar.update(min(estimated_stars, remaining))
                    pbar.set_postfix({
                        'skipped': skipped, 
                        'downloaded': downloaded, 
                        'failed': failed,
                        'stars': f'{total_stars_estimated:,}'
                    })
                else:
                    file_size_mb = filepath.stat().st_size / (1024 * 1024)
                    print(f"  ✓ Downloaded {filename} ({file_size_mb:.1f} MB, ~{estimated_stars:,} stars)")
                    sys.stdout.flush()
            else:
                failed += 1
                if HAS_TQDM:
                    pbar.set_postfix({
                        'skipped': skipped, 
                        'downloaded': downloaded, 
                        'failed': failed,
                        'stars': f'{total_stars_estimated:,}'
                    })
                else:
                    print(f"  ✗ Failed: {filename}")
                    sys.stdout.flush()
            
            # Small delay to be nice to the server
            time.sleep(0.1)
            
    except KeyboardInterrupt:
        print("\n\nDownload interrupted by user.")
        if HAS_TQDM and pbar:
            pbar.close()
        raise
    finally:
        if HAS_TQDM and pbar:
            pbar.close()
    
    print(f"\nDownload complete!")
    print(f"  Downloaded: {downloaded} files")
    print(f"  Skipped: {skipped} files")
    print(f"  Failed: {failed} files")
    print(f"  Estimated stars: {total_stars_estimated:,} (target: {star_limit:,})")

if __name__ == "__main__":
    main()
