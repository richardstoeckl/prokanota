from pathlib import Path


def _read_version() -> str:
    version_file = Path(__file__).with_name("prokanota.VERSION")
    return version_file.read_text(encoding="utf-8").strip()


__version__ = _read_version()
