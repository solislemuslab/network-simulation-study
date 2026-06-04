#!/bin/bash

# DAGMan passes the dataset ID
dataset=$1

# CD into the specific job folder
JOB_DIR="/home/justison/jobs/job_${dataset}"
cd $JOB_DIR || exit 1



# SAFETY CHECK: Only delete the raw files if the tarball was successfully created!
if [ -f "emp_est_out${dataset}.tar.gz" ]; then
    find . -maxdepth 1 -type f -not -name '*.tar.gz' -delete
else
    exit 1
fi