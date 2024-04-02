#!/usr/bin/env python

# Copyright 2008-2024 Michael M. Hoffman <michael.hoffman@utoronto.ca>

from setuptools import Extension, setup
from pathlib import Path

AWFMINDEX_BUILD_DIR = (Path("AvxWindowFmIndex") / "build").resolve()

generate_index_source_files = ["src/newmap-generate-index.c"]
count_kmers_source_files = ["src/newmap-count.c"]

dynamic_libraries = ["gomp"]  # for OpenMP

# Libraries exist in git submodules
library_dirnames = list(map(str, [AWFMINDEX_BUILD_DIR]))
include_dirnames = list(map(str, [AWFMINDEX_BUILD_DIR]))
static_libraries = list(map(str,
                        [AWFMINDEX_BUILD_DIR / "libawfmindex_static.a",
                         AWFMINDEX_BUILD_DIR / "libfastavector_static.a",
                         AWFMINDEX_BUILD_DIR / "libdivsufsort64.a"]))

# Build against the stable ABI to support 3.9 onwards
c_define_macros = [("Py_LIMITED_API", "0x03090000")]

generate_index_module = Extension(
    '_c_newmap_generate_index',  # needs to match C file PyInit definition
    sources=generate_index_source_files,
    include_dirs=include_dirnames,
    library_dirs=library_dirnames,
    libraries=dynamic_libraries,
    extra_objects=static_libraries,
    define_macros=c_define_macros,
    py_limited_api=True)

count_kmers_module = Extension(
    '_c_newmap_count_kmers',  # needs to match C file PyInit definition
    sources=count_kmers_source_files,
    include_dirs=include_dirnames,
    library_dirs=library_dirnames,
    libraries=dynamic_libraries,
    extra_objects=static_libraries,
    define_macros=c_define_macros,
    py_limited_api=True)


if __name__ == "__main__":
    # place extension in the base umap package
    setup(ext_package="newmap",
          ext_modules=[generate_index_module,
                       count_kmers_module])
