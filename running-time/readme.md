# SNaQ running time study

Goal: Study the variables that play a role in the running time of SNaQ.

## Data
- We will use the simulated data for snaq and ransanec projects (`knownGT`)
- We have 5 networks of increasing number of tips: n6, n10, n15, n25, n50
- For each network, we have 10 replicates obtained as follow:
    - 2 replicates per number of gene trees: 30,100,300,1000,3000 (`1_astral.in` and `2_astral.in` renamed as `1_astral_30.in` for 30 gene trees e.g.)
    - we do not expect to see changes in running time due to number of gene trees

For n6, n10, n15:
```shell
cd Dropbox/Documents/solislemus-lab/lab-members/grad-student-projects/rotations/nathan-projects/ransanec/ransanec/data/simulations/knownGT

cp n15_n3000/1_astral.in /Users/Clauberry/Dropbox/Documents/solislemus-lab/lab-members/postdoc-projects/gustavo-projects/network-simulation-study/running-time/data/n15

mv /Users/Clauberry/Dropbox/Documents/solislemus-lab/lab-members/postdoc-projects/gustavo-projects/network-simulation-study/running-time/data/n15/1_astral.in /Users/Clauberry/Dropbox/Documents/solislemus-lab/lab-members/postdoc-projects/gustavo-projects/network-simulation-study/running-time/data/n15/1_astral_3000.in 

cp n15_n3000/2_astral.in /Users/Clauberry/Dropbox/Documents/solislemus-lab/lab-members/postdoc-projects/gustavo-projects/network-simulation-study/running-time/data/n15

mv /Users/Clauberry/Dropbox/Documents/solislemus-lab/lab-members/postdoc-projects/gustavo-projects/network-simulation-study/running-time/data/n15/2_astral.in /Users/Clauberry/Dropbox/Documents/solislemus-lab/lab-members/postdoc-projects/gustavo-projects/network-simulation-study/running-time/data/n15/2_astral_3000.in 
```

## Running SNaQ

(maybe we want to run in on chtc)

1. Create CF tables
2. Run snaq in parallel for h=0 and nruns=1