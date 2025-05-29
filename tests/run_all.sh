#!/usr/bin/env bash

newmap index \
    --compression-ratio=32 \
    --seed-length=1 \
    data/genome.fa

newmap search \
    --output-directory=unique_lengths \
    --search-range=4:10 \
    data/genome.fa

newmap track unique_lengths/*
newmap track --multi-read - 8 unique_lengths/*
