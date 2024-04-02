#!/usr/bin/env bash
# This script builds the AvxWindowFmIndex library and it's dependencies

# Fetch all submodules if they don't exist
git submodule update --init --recursive

# Build AvxWindowFmIndex libraries
cd AvxWindowFmIndex

# Force no shared libraries
# Force position independent code for creating shared library

# Debug build
# cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON .

# Release build
cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON .
make
