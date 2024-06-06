import gzip
from pathlib import Path
import sys


def optional_gzip_open(file_path: Path, mode: str):
    """If the filename ends with .gz, return a gzip file object,
    otherwise return a regular file object."""

    if file_path.suffix == ".gz":
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def verbose_print(verbose: bool, *args):
    if verbose:
        print(*args, file=sys.stderr)
