#!/bin/bash

# run Julia script

tar -xzvf files.tar.gz

julia htc_network_estimation.jl ./ ./out/ 1

tar -czf output_$1.tar.gz out/*

