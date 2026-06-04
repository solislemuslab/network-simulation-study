#!/bin/bash
N_CPUS=$1
JULIA_SCRIPT=$2


# Unzip your raw data (assuming files.tar.gz contains the raw inputs)
tar -xzvf "files.tar.gz"


# Run the CF simulation and tree generation
julia $JULIA_SCRIPT $N_CPUS
