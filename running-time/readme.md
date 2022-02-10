# SNaQ running time study

Goal: Study the variables that play a role in the running time of SNaQ.

## Data
- We will use the simulated data for snaq and ransanec projects (`knownGT`)
- We have 5 networks of increasing number of tips: n6, n10, n15, n25, n50
- For each network (n6,n10,n15), we have 10 replicates obtained as follow:
    - 2 replicates per number of gene trees: 30,100,300,1000,3000 (`1_astral.in` and `2_astral.in` renamed as `1_astral_30.in` for 30 gene trees e.g.)
    - we do not expect to see changes in running time due to number of gene trees
- For each network (n25, n50), we have 10 replicates obtained as follow:
    - 2 replicates per number of gene trees: 50, 100, 500, 1000, 5000 e.g. `n25h5-gt50-1.tre`, `n25h5-gt50-2.tre`
    - Note that files with 5000 gene trees had to be saved into 5 files of 1000 gene trees.


## Running SNaQ

(maybe we want to run it on chtc)

1. Create CF tables (save running time too and memory also)
2. Run snaq in parallel for h=0 and nruns=1
