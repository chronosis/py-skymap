# Py-Skymap

Generate galactic hemisphere maps centered on any target star using Gaia DR3 data. This tool creates beautiful polar projection maps showing the night sky as it would appear from the perspective of distant stars, including known bright stars, galaxies, and deep-sky objects.

## Features

- **Gaia DR3 Integration**: Downloads and caches star data from the European Space Agency's Gaia mission
- **SQLite Caching**: Efficient local caching of star data to avoid repeated API calls
- **Custom Star Perspectives**: View the sky from any star's perspective, not just Earth
- **Known Objects**: Includes 100+ bright stars, galaxies (LMC, SMC, Andromeda, Triangulum), and deep-sky objects (nebulas, clusters)
- **Dynamic Reference**: Automatically includes Sol (the Sun) as a reference point when viewing from other stars
- **High-Quality Output**: Generates 24x24 inch PNG images with proper star magnitudes and object labeling

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
poetry run python skymap-gen.py Aldeberan

# Generate maps with 100,000 stars in cache
poetry run python skymap-gen.py Aldeberan 100000

# Generate maps for Sol (the Sun)
poetry run python skymap-gen.py Sol

# Generate maps for Proxima Centauri
poetry run python skymap-gen.py "Proxima Centauri"
```

### Command-Line Options

- `target_star`: Name of the target star (required). Can be:
  - HD catalog number (e.g., `HD100000`)
  - Common name (e.g., `Sol`, `Proxima Centauri`, `Alpha Centauri`)
  - Any name resolvable by Gaia Archive
  
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

### Output

The script generates two PNG files in the `./images/` directory:
- `images/<target_star>_north_hemisphere.png` - Northern galactic hemisphere
- `images/<target_star>_south_hemisphere.png` - Southern galactic hemisphere

The `images/` directory is created automatically if it doesn't exist. Images are 24x24 inches at 150 DPI, suitable for printing or detailed viewing.

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

## Project Structure

```
py-skymap/
├── skymap-gen.py          # Main script for generating sky maps
├── download-gaia-csv.py   # Script for downloading CSV files from ESA CDN
├── gaia-query.py          # Utility script for Gaia queries
├── pyproject.toml         # Poetry dependencies configuration
├── .gitattributes         # Git LFS configuration for PNG files
├── gaia_cache/            # SQLite database cache (git-ignored)
│   └── gaia_cache.db      # Cached star data
└── images/                # Generated sky map images (tracked with Git LFS)
    └── *.png
```

## Dependencies

Managed via Poetry (see `pyproject.toml`):

- **numpy**: Numerical computations
- **matplotlib**: Plotting and visualization
- **astropy**: Astronomical coordinate transformations and calculations
- **astroquery**: Querying Gaia Archive
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
