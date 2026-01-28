try:
    from tqdm import tqdm as _tqdm  # type: ignore[assignment]
    HAS_TQDM: bool = True
except ImportError:  # pragma: no cover - fallback path
    HAS_TQDM = False

    def _tqdm(iterable, **kwargs):
        return iterable


# Public export used by callers
tqdm = _tqdm

