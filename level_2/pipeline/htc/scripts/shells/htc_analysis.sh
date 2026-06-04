#!/bin/bash

# run Julia script
mkdir out

for f in *.tar.gz
do
  tar -xzvf "$f"
done


julia analysis_single.jl

tar -czf $2_out${1}.tar.gz clusters.csv found_clades.csv squirrel.csv
