import gzip
from pathlib import Path
from typing import Any, IO, TextIO, Union
import sys


def optional_gzip_open(file_path: Path,
                       mode: str) -> Union[gzip.GzipFile, TextIO, IO[Any]]:
    """If the filename ends with .gz, return a gzip file object,
    otherwise return a regular file object."""

    if file_path.suffix == ".gz":
        return gzip.open(file_path, mode)  # GzipFile | TextIO
    else:
        return open(file_path, mode)  # IO[Any]


def verbose_print(verbose: bool, *args):
    if verbose:
        print(*args, file=sys.stderr)
