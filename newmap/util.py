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

# NB: Must use ceil when attempting to go to higher lengths
# During binary search
def ceil_div(a, b):
    """Returns the integer ceiling of a / b."""
    # NB: math.ceil(a / b) is not the same as this function
    # And is floating point-based
    return -(a // -b)
