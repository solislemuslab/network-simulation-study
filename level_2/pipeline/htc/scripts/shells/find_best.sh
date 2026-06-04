#!/bin/bash

# Arguments passed by DAGMan
h_val=$1
dataset=$2

# DAGMan runs this from the submit directory, so we must CD into the specific job folder
JOB_DIR="/home/justison/jobs/job_${dataset}"
cd $JOB_DIR || exit 1

echo "Finding best network for dataset ${dataset} at hmax=${h_val}..."

# ------------------------------------------------------------------
# 1. grep: Strictly finds the line containing "-Ploglik ="
# 2. awk: Splits the filename off the front, and grabs the last word (the score)
# 3. sort: Sorts the scores numerically (lowest score floats to the top)
# 4. head: Grabs the winning top line
# 5. awk: Extracts just the filename of the winner
# ------------------------------------------------------------------
best_file=$(grep -H "\-Ploglik =" h_${h_val}_run_*.out | awk '{split($1,a,":"); print a[1], $NF}' | sort -k2 -n | head -n 1 | awk '{print $1}')

# Safety check: Ensure we actually found a file
if [ -z "$best_file" ]; then
    echo "ERROR: Could not find any valid SNaQ .out files for h=${h_val}!"
    exit 1
fi

echo "Winning file is: $best_file"

# Secure the winner to the exact name the next DAG step is looking for
cp "$best_file" "h_${h_val}.out"

# ------------------------------------------------------------------
# THE CLEANUP
# ------------------------------------------------------------------
rm h_${h_val}_run_*.out

echo "Cleanup complete. Best network secured for h=${h_val}."