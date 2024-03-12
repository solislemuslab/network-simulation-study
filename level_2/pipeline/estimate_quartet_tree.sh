#!/usr/bin/bash

## positional unnamed arguments
GENETREES=$1
OUTFILE=$2

treeqmc --fast -i $GENETREES -o $OUTFILE
