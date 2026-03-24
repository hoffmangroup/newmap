#!/usr/bin/env bash

set -x

newmap --version

newmap index \
    --compression-ratio=32 \
    --seed-length=1 \
    data/genome.fa

newmap search \
    --verbose \
    --output-directory=unique_lengths \
    --search-range=4:10 \
    data/genome.fa

newmap track --verbose unique_lengths/*
newmap track --verbose --multi-read - 8 unique_lengths/*
