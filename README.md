# Py-Skymap

Generate galactic hemisphere maps centered on any target star using Gaia DR3 data. This tool creates beautiful polar projection maps showing the night sky as it would appear from the perspective of distant stars, including known bright stars, galaxies, and deep-sky objects.

## Features

- **Gaia DR3 Integration**: Downloads and caches star data from the European Space Agency's Gaia mission
- **SQLite Caching**: Efficient local caching of star data to avoid repeated API calls
- **Custom Star Perspectives**: View the sky from any star's perspective, not just Earth
- **Known Objects**: Includes 100+ bright stars, galaxies (LMC, SMC, Andromeda, Triangulum), and deep-sky objects (nebulas, clusters)
- **Dynamic Reference**: Automatically includes Sol (the Sun) as a reference point when viewing from other stars
- **High-Quality Output**: Generates PNG images at 24"×24" @ 150 DPI by default, with proper star magnitudes and object labeling
- **Proportional Labeling**: Text sizes scale with image dimensions; semi-transparent black stroke for legibility on varied backgrounds
- **Stellar Color Mapping**: Star colors are derived from Gaia's BP-RP color index, reflecting each star's actual stellar class (effective temperature)—blue for hotter stars, red for cooler
- **Depth-Ordered Rendering**: Stars are drawn back to front (farthest first, nearest last) so nearer stars correctly occlude more distant ones
- **Magnitude-Based Sizing**: Point sizes reflect each star's apparent magnitude from the target star's perspective—brighter stars appear larger

## Requirements

- Python 3.11 or higher
- Poetry (for dependency management)
- Git LFS (for tracking large PNG files)

## Installation

### 1. Install Poetry

If you don't have Poetry installed, follow the [official installation guide](https://python-poetry.org/docs/#installation).

### 2. Install Git LFS

```bash
# On Ubuntu/Debian
sudo apt-get install git-lfs

# On macOS (with Homebrew)
brew install git-lfs

# On Windows (with Chocolatey)
choco install git-lfs
```

### 3. Clone and Setup

```bash
# Clone the repository
git clone <repository-url>
cd py-skymap

# Initialize Git LFS (if not already done)
git lfs install

# Install dependencies with Poetry
poetry install
```

## Usage

### Basic Usage

Generate sky maps for a target star:

```bash
poetry run python skymap-gen.py <target_star> [stars]
```

**Examples:**

```bash
# Generate maps with default 50,000 stars in cache
poetry run python skymap-gen.py Aldebaran

# Generate maps with 100,000 stars in cache
poetry run python skymap-gen.py Aldebaran 100000

# Generate maps for Sol (the Sun)
poetry run python skymap-gen.py Sol

# Generate maps for Proxima Centauri (quote names with spaces)
poetry run python skymap-gen.py "Proxima Centauri"

# Generate maps for multiple stars (comma-separated)
poetry run python skymap-gen.py "Sol,Aldebaran,Proxima Centauri" 100000
```

### Command-Line Options

- `target_star`: Name of the target star(s) (required). Can be:
  - HD catalog number (e.g., `HD100000`)
  - Common name (e.g., `Sol`, `Proxima Centauri`, `Alpha Centauri`)
  - Any name resolvable by Gaia Archive
  - **Multiple stars**: Comma-separated list; wrap the whole list in quotes. Quote individual names with spaces (e.g., `"Sol,Aldebaran,Proxima Centauri"`)
  
- `stars`: Number of stars to download from Gaia (optional, default: 50,000)
  - Higher numbers provide more detail but take longer to download and process
  - Recommended: 50,000-200,000 for most use cases

- `--force-refresh`: Force re-download of cached data from Gaia
  ```bash
  poetry run python skymap-gen.py HD118246 --force-refresh
  ```

- `--dump-positions`: Write 3D star positions to database without generating images
  ```bash
  poetry run python skymap-gen.py Sol 50000 --dump-positions
  ```

- `--magnitude-limit`: Faintest magnitude to render (overrides default from constants)
  ```bash
  poetry run python skymap-gen.py Sol 50000 --magnitude-limit 13
  ```

- `--point-size-min`: Smallest star point size in the plot (default: 3)
  ```bash
  poetry run python skymap-gen.py Sol 50000 --point-size-min 2
  ```

- `--point-size-max`: Largest star point size in the plot (default: 132)
  ```bash
  poetry run python skymap-gen.py Sol 50000 --point-size-max 200
  ```

### Output

The script generates three PNG files per target star in the `./images/` directory (spaces in the target name become dashes in filenames):
- `images/<target_star>_north_hemisphere.png` - Northern galactic hemisphere
- `images/<target_star>_south_hemisphere.png` - Southern galactic hemisphere
- `images/<target_star>_360_degree.png` - Full 360° view (east and west hemispheres combined)

The `images/` directory is created automatically if it doesn't exist. By default, images are 24"×24" at 150 DPI (3600×3600 pixels), suitable for printing or detailed viewing.

### Example Images

**Sol (Sun) — Earth's perspective:**

| North hemisphere | South hemisphere |
|------------------|------------------|
| [![Sol North](images/SOL_north_hemisphere.png)](images/SOL_north_hemisphere.png) | [![Sol South](images/SOL_south_hemisphere.png)](images/SOL_south_hemisphere.png) |

**Aldebaran — from ~65 light-years away:**

| North hemisphere | South hemisphere |
|------------------|------------------|
| [![Aldebaran North](images/Aldebaran_north_hemisphere.png)](images/Aldebaran_north_hemisphere.png) | [![Aldebaran South](images/Aldebaran_south_hemisphere.png)](images/Aldebaran_south_hemisphere.png) |

**More examples** in the `images/` folder (each target includes north, south, and 360° views):

- **[Alpha Centauri](images/Alpha-Centauri_north_hemisphere.png)** — Relatively nearby star to Sol (part of the closest stellar system to the Sun).
- **[Betelgeuse](images/Betelgeuse_north_hemisphere.png)** — Red supergiant over 100 ly from Earth (in Orion).
- **[HD100000](images/HD100000_north_hemisphere.png)** (BD-22 3152) — Giant star ~1000 ly from Earth, closer to the galactic center.
- **[HD118246](images/HD118246_north_hemisphere.png)** — Star high above the main galactic plane.
- **[Polaris](images/Polaris_north_hemisphere.png)** — The North Star; appears almost directly above Earth's north pole.
- **[Proxima Centauri](images/Proxima-Centauri_north_hemisphere.png)** — Relatively nearby star to Sol (closest known star to the Sun).

*These example images were generated using ~101M stars (about 10% of Gaia DR3) in the local cache (magnitude ≤16) and a visibility limit of 11.5. The initial download of 101M stars from Gaia took approximately 10 hours; generating the skymaps from cached data takes about 5 minutes per target on a decently powered system.*

## Data Management

### SQLite Cache

Star data is cached locally in `gaia_cache/gaia_cache.db` to avoid repeated API calls. The cache includes:
- Star positions (RA, Dec, parallax)
- Magnitudes (photometric bands)
- Source IDs for deduplication

### Cache Management

- First run: Downloads data from Gaia Archive (may take time depending on star limit)
- Subsequent runs: Uses cached data (much faster)
- Force refresh: Use `--force-refresh` to re-download data

### Downloading CSV Files (Alternative)

For bulk data downloads, use the `download-gaia-csv.py` script:

```bash
poetry run python download-gaia-csv.py [stars]
```

This downloads CSV.gz files directly from the ESA CDN, bypassing the Gaia API.

### Bulk Download from VizieR

You can bulk-download star data from [VizieR](https://vizier.cds.unistra.fr/) into the same SQLite cache. Supported for **Gaia DR3** (catalog `I/355/gaiadr3`): rows are merged into `gaia_source` so `skymap-gen` can use them like the Gaia Archive or CSV pipeline.

```bash
# Default: Gaia DR3, 50,000 rows
poetry run python download-vizier.py

# Gaia DR3 with a custom row limit
poetry run python download-vizier.py I/355/gaiadr3 100000
poetry run python download-vizier.py I/355/gaiadr3 --limit 50000

# Full catalog (very large for Gaia DR3)
poetry run python download-vizier.py I/355/gaiadr3 --all
```

Use `--no-merge` to fetch without writing to the cache (e.g. for testing). Other VizieR catalogues are not yet mapped into `gaia_source`; only Gaia DR3 is supported for merging.

### Simbad Cache (filling gaps)

Simbad data can be queried and cached in the same SQLite database (`gaia_cache.db`) to flesh out gaps—e.g. objects not in Gaia, or updated positions/magnitudes/parallax for named objects. The cache table `simbad_cache` stores `main_id`, `ra`, `dec`, `parallax_mas`, `vmag`, `otype`, and `cached_at`.

- **Populate the cache** using the bright-star list (or your own list):
  ```bash
  poetry run python fetch-simbad-cache.py
  poetry run python fetch-simbad-cache.py "M1,M2,NGC 1976"
  poetry run python fetch-simbad-cache.py --names "Sirius,Aldebaran,Vega"
  ```
- **API**: Use `lib.simbad_client.query_simbad_and_cache(cache_db, identifiers)` to fetch and cache in bulk, or `get_from_simbad_or_cache(cache_db, identifier)` for a single lookup (cache-first). Read back with `lib.sqlite_helper.get_simbad_from_cache(db_path, identifier)` or `get_all_simbad_cached(db_path)`.

Simbad is rate-limited (recommended 5–10 queries/sec); the client batches requests and adds a short delay between them.

**Bulk Simbad via TAP:** Simbad does not provide a full catalog dump like Gaia’s CDN; it’s a dynamic database. For bulk data you can use **TAP (Table Access Protocol)** with ADQL:

- **astroquery**: `Simbad.query_tap("SELECT TOP 50000 main_id, ra, dec FROM basic WHERE ...", async_job=True, timeout=300, maxrec=2000000)` for large queries. Use `async_job=True` and a `timeout` to avoid timeouts.
- **Cone search**: restrict by sky region with `CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', ra_center, dec_center, radius_deg)) = 1`.
- **Criteria**: filter by `otype`, or join with the `flux` table for magnitude limits. See [SIMBAD ADQL examples](https://cds.unistra.fr/help/documentation/simbad-more/adql-simbad) and the [ADQL cheat sheet](https://simbad.u-strasbg.fr/simbad/tap/help/adqlHelp.html).

You can run such a TAP query, then insert the result rows into `simbad_cache` (with the same columns: main_id, ra, dec, parallax_mas, vmag, otype) to bulk-fill the cache. For complete published catalogues, CDS recommends [VizieR](https://vizier.cds.unistra.fr/) rather than Simbad.

## Project Structure

```
py-skymap/
├── skymap-gen.py          # Main script for generating sky maps
├── fetch-simbad-cache.py  # Populate Simbad cache (optional, for filling gaps)
├── download-gaia-csv.py   # Script for downloading CSV files from ESA CDN
├── download-vizier.py     # Bulk download from VizieR (e.g. Gaia DR3 → gaia_source)
├── gaia-query.py          # Utility script for Gaia queries
├── pyproject.toml         # Poetry dependencies configuration
├── .gitattributes         # Git LFS configuration for PNG files
├── gaia_cache/            # SQLite database cache (git-ignored)
│   └── gaia_cache.db      # Cached Gaia + Simbad + VizieR data
├── lib/                   # Shared library (gaia_client, simbad_client, vizier_client, etc.)
└── images/                # Generated sky map images (tracked with Git LFS)
    └── *.png
```

## Dependencies

Managed via Poetry (see `pyproject.toml`):

- **numpy**: Numerical computations
- **matplotlib**: Plotting and visualization
- **astropy**: Astronomical coordinate transformations and calculations
- **astroquery**: Querying Gaia Archive and Simbad
- **tqdm**: Progress bars
- **requests**: HTTP requests
- **beautifulsoup4**: HTML parsing

## How It Works

1. **Target Star Resolution**: Looks up the target star in the cache or queries Gaia Archive
2. **Data Download**: Downloads star data from Gaia DR3 (or uses cache)
3. **Coordinate Transformation**: Transforms star positions from Earth-centric to target-star-centric coordinates
4. **Magnitude Calculation**: Calculates apparent magnitudes from the target star's perspective
5. **Hemisphere Splitting**: Separates stars into northern and southern hemispheres based on z-component
6. **Known Objects**: Adds bright stars, galaxies, and deep-sky objects with proper distances
7. **Visualization**: Creates polar projection maps with proper scaling and labeling

## Known Objects Included

### Galaxies
- Large Magellanic Cloud (LMC)
- Small Magellanic Cloud (SMC)
- Andromeda Galaxy (M31)
- Triangulum Galaxy (M33)

### Bright Stars (100+)
Includes well-known stars like:
- Rigel, Betelgeuse, Antares, Aldebaran
- Pollux, Polaris, Sirius, Vega
- And many more...

### Deep-Sky Objects
- Nebulas (Orion, Carina, etc.)
- Globular clusters
- Open clusters (Pleiades, Hyades, etc.)

All objects are positioned using their actual distances/parallax values when available.

## Troubleshooting

### "ModuleNotFoundError: No module named 'numpy'"

Ensure dependencies are installed:
```bash
poetry install
```

### "Job timeout/aborted" from Gaia Archive

- Reduce the number of stars (use a lower value for the `stars` parameter)
- Try again later (Gaia Archive may be experiencing high load)

### Images appear empty or sparse

- Increase the number of stars (try 100,000 or 200,000)
- Check that the cache has data: `ls -lh gaia_cache/gaia_cache.db`

### Git LFS issues

If PNG files aren't tracked properly:
```bash
git lfs install
git add .gitattributes
git add images/*.png
```

## License

See [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Gaia Mission**: European Space Agency's Gaia space observatory
- **Astropy**: Astronomical Python library
- **Astroquery**: Astronomical data querying tools

## Contributing

Contributions are welcome! Please ensure:
- Code follows Python best practices
- Dependencies are added to `pyproject.toml`
- Large files (PNGs) are tracked with Git LFS
- Tests pass (if applicable)

## References

- [Gaia Archive](https://gea.esac.esa.int/archive/)
- [Gaia DR3 Documentation](https://www.cosmos.esa.int/web/gaia/dr3)
- [Astropy Documentation](https://docs.astropy.org/)
