This folder contains the scripts used to simulate data, estimate networks, and summarize/compare the results.



## Data Generation

There are two main simulation sets that are run:

### Main Analysis
This consists of 36 different parameter sets of 150 simulated networks and their concordance factors.
The different parameterizations are as follows
- *ntaxa* - 15,20 ,30
- *ngenes* - 100, 1000, 10000
- *hybridization rate* - 0.02, 0.04
- *restricted to level-1* - yes, no

There are three main levels to the data generation of this analysis:
- par_no - this is the parameter ID that corresponds to a specific parameter setting
- phy_no - for each parameterization, 150 networks were generated, this is the network ID
- rep_no - for each network, 30 replicate sets of CFs were generated, this corresponds to the CFs ID

to simplify running analyses on HTC, these three different levels get flattened into a single job ID. `01.generate_datasets.R` creates a file called `job_map.csv` that maps job IDS to the respective par_no, phy_no, and rep_no values

The following scripts in the `pipeline` folder are used to generate these data:
- `00.generate_seeds.R` - This file contains the main simulation/inference settings, sets the global seed, creates seeds at varies points in the data simulation/inference, and creates the `pars.csv` file containing information that maps parameterizations to a specific parameter ID. 

- `01.generate_datasets.R` - This is the main file for data simulation and prep. It first generates phylogenetic networks under a BDH process. It will run until generating the specified number of networks that meet the criteria: Are level-1 (or at least level-2 depending on the parameterization) and have no 2-cycles, 3-cycles, and the root is the last stable ancestor. This script will call `02.sim_gts_calc_cf.jl` (see below), save relevant data in the `data` folder and then compress all necessary data into the `jobs` folder to be sent to HTC. 

- `02.sim_gts_calc_cf.jl` This file reads the individual networks, simulates gene trees under the NMSC, computes concordance factors and uses tree-QMC to estimate a starting tree for SNaQ. For this file to run the user must have tree-QMC available and it must be callable from the command line as `tree-qmc`

At the end of this process, the `jobs` folder will be populated with `job_i` directories. Each of these folders will contain a `files.tar.gz` tarball. In this compressed folder there are 3 files:
- `CFs.csv`: the concordance factors
- `starting_tree.newick`: the starting tree for SNaQ
- `est_pars.csv`: a CSV containing the relevant parameter values and seeds for SNaQ.

Some information is also saved in the data folder. The directory structure uses the par_no, phy_no, and rep_no structure to store appropriate information at each level. For example each simulated phylogenetic network will be stored at `./data/pars_i/net_j/network.extnewick` while CFs are different at each rep_no and will be stored as `./data/pars_i/net_j/rep_k/CFs.csv`

The `jobs` directory is sent to HTC to be analyzed (see below)

### Toy Analysis
This analysis uses some of the networks generated under the **Main Analysis** but sets their inheritance proportions to 0.5 before simulating data. It assumes phylogenetic networks were already generated with `01.generate_datasets.R` and uses modified scripts to generate the CFS and prep jobs for the HTC.
- `generate_toy_datasets.R` This is similar to the `01.generate_datasets.R` file but doesn't actually simulate networks. It calls `02.sim_gts_toy.jl` and then packages all starting data into a `jobs_toy` folder for the HTC
- `02.sim_gts_toy.jl` to generate CFS and starting trees. 


## HTC Network inference and Analysis 

HTC resources were used for the main SNaQ inference and other computationally intense analyses. See the README in the `htc_files` folder for a description of the files for inference and analysis. 


## Network Inference and Comparison

Here we analyze the properties of the simulated networks, estimated networks, and compare these results. The following files are used here summarize information on the simulated networks:

`sim_sums.R`: This computes some basic information about the network such as the level, number of blobs, class membership (tree-child, tree-based, FU-stable), and confirms the lack of small cycles. Results here are saved to `./data/sim_net_props.csv`

`sim_sums.jl`: This script computes more information about the network, computes specific properties about individual hybridizations, and branch lengths. it will generate the following files and folders:
- - `data/pars_i/hyb_dat.csv`: This contains information about each hybridization found across the networks within the parameterization
- - `data/pars_i/branch_lengths.csv`: This contains the internal branch lengths found across the networks 
- - `data/pars_i/net_j/hybs/H_k/stacked.csv`: This contains the hybrid nodes that are stacked above or below the hybridization
- - `data/pars_i/net_j/hybs/H_k/cycle_nodes.csv`: The nodes that are part of the cycle of the hybridization
- - `data/pars_i/net_j/hybs/H_k/cycle_hybrids.csv`: the nodes that are a part of the cycle of the hybridization and denotes whether they are above or below on the cycle

`sim_sums_reps`: This script breaks a network down into its subnetworks and computes the likelihood of the data given each subnetwork. This script is run at the rep_no level since it uses CFs. NOTE that this script is slow to run and has been relegated to being performed on the HPC. Assuming things ran well, it generates a csv of the likelihood at each level-1 subnetwork in `./data/pars_i/phy_j/rep_k/liks.csv`.

Note that SNaQ inference depends on the maximum number of hybridizations (hmax) allowed during the search. We ran SNaQ for values of hmax 0-5. This adds another level to our nested dataset, h_no. Generally information is summarizedacross par_no settings (i.e., the csv will be in the par_i folder but information on phy_no, rep_no, and h_no are contained within the csv). it generates the following files in the newly created `summarized_results`.

These files summarize the estimated networks to the true network:

- `net_properties.R`: This computes the tree of blobs for both the estimated and true networks and compares them. It creates the following files:
- - `./summarized_results/par_no/tob.csv` This computes distances (various metrics) between the estimated and true tree of blobs
- - `./summarized_results/par_no/tob/compat.csv` This contains information on whether the individual blobs found within the estimated network are compatible with any of the blobs in the true network
- `est_net_gof.jl`: This performs a goodness of fit test for all of the estimated networks. Results are saved in `./summarized_results/par_i/GoF.csv`.
- `analysis.jl`:This script finds the number of 'found' hybrids and clusters in the network, in both compatible and strict senses. It also computes the HWCD between the estimated net and the true network and its subnetworks. It generates three files:
- - `../summarized_results/pars_i/clusters.csv`: this contains summary information about how many hybrids/hybrid clusters were found, HWCD. 
- - `../summarized_results/pars_i/found_clades.csv`: `This contains specific information about *which* clusters were found in the estimated network
- - `../summarized_results/pars_$par_no/subnet_dat.csv`: This file contains information about the distances between the estimated network and the true subnetworks and their likelihood differences


These files summarize process and summarize the information into figures and linear models
`process_output.R`: this takes the data from the HPC and puts it in a more workable folder structure.

`results_summary.R`: filters the data based on a few different criteria and then feeds the filtered datasets to `filtered_sum.R` to make figures and compute summary statistics





