#!/bin/bash

echo "--- DEBUG INFO ---"
echo "I am running on host: $(hostname)"
echo "My arguments are: $@"
echo "Arg 1 (CPUs): $1"
echo "Arg 2 (h_val): $2"
echo "Arg 3 (job_id): $3"
echo "Arg 4 (start tree): $4"
echo "Arg 5 (prefix): $5"
echo "Arg 6 (run_no): $6"
echo "------------------"


# run Julia script
mkdir out

for f in *.tar.gz
do
  tar -xzvf "$f"
done

julia emp_est_single.jl $1 $2 $3 $4 $5 $6


#tar -czf $5_out${3}.tar.gz  -C out .

