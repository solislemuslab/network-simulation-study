This folder contains the scripts used to simulate data, estimate networks, and summarize/compare the results.

`aux_functions.jl` and `functions.R` contain helper functions for computing stats and bookkeeping.

The analysis files use 


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
- par_no - this is the parameter ID that corresponds to a specific parameter setting. The specific parameter setting can be found in `/data/pars.csv`
- phy_no - for each parameterization, 150 networks were generated, this is the network ID
- rep_no - for each network, 30 replicate sets of CFs were generated, this corresponds to the CFs ID

to simplify running analyses on HTC, these three different levels get flattened into a single job ID.
The `get_jobID` function defined in the `functions.R` define the mapping from parameter numbers parameter number to flat job ID 
 `01.generate_datasets.R` creates a file called `job` that maps job IDS to the respective par_no, phy_no, and rep_no values. The function `get_jobID` from `functions.R` will also achieve this goal.

The following scripts in the `pipeline` folder are used to generate these data:
- `00.generate_seeds.R` - This file contains the main simulation/inference settings, sets the global seed, creates seeds at varies points in the data simulation/inference, and creates the `pars.csv` file containing information that maps parameterizations to a specific parameter ID. 

- `01.generate_datasets.R` - This is the main file for data simulation and prep. It first generates phylogenetic networks under a BDH process. It will run until generating the specified number of networks that meet the criteria: Are level-1 (or at least level-2 depending on the parameterization) and have no 2-cycles, 3-cycles, and the root is the last stable ancestor. This script will then compress all necessary data into a tar.gz file called `job_ID` that contains all nessecary downstream material to be used in simulation, inference and analysis. 

At the end of this process, the `jobs` folder will be populated with `job_i` directories. Each of these folders will contain a `files.tar.gz` tarball. In this compressed folder there are 3 files:
- `starting_tree.newick`: the starting tree for SNaQ
- `est_pars.csv`: a CSV containing the relevant parameter values and seeds for SNaQ.

Some information is also saved in the data folder. The directory structure uses the par_no, phy_no, and rep_no structure to store appropriate information at each level. For example each simulated phylogenetic network will be stored at `./data/pars_i/net_j/network.extnewick` while CFs are different at each rep_no and will be stored as `./data/pars_i/net_j/rep_k/CFs.csv`

The `jobs` directory is sent to HTC to be analyzed (see below)

The `dag_maker.jl` script scans output directories to detect missing datasets across various settings, phylogenies, and replicates and builds an HTCondor Directed Acyclic Graph (DAG) workflow file (`master_workflow.dag`).
It packages missing job directories into a compressed archive (`subset_jobs.tar.gz`). This archive was sent to htc resources for processing


## HTC Network inference and Analysis 

HTC resources were used for the main SNaQ inference and other computationally intense analyses.
See the README in the `htc_files` folder for a description of the files for inference and analysis. 

These analyses use modified versions of `PhyloNetworks.jl` and `SNaQ.jl`.
The modified version of `SNaQ` is avaliable in `julia_est.sif` which contains a development version of SNaQ 2.0 that was used to speed up analysis times.
Specifically, the sif is built from commit 65d4165d802609afb05e49305c3249d093c3843b on the treechild-galled branch.
The modified version of `PhyloNetworks` is in `julia_analysis.sif` and contains helper functions for computing some summaries of the networks.
The modified version primarilty has additional functions defined the in the file `network_reductions.jl`

## Network and Comparison

The analyses in Julia use a modified version of `PhyloNetworks.jl` with some additional functionality to compute the summary statistics. 
This version is avaliable directly in the `julia_analysis.sif` container. 

`process_output.R`: this takes the data from the HPC and puts it in a more workable folder structure.

NOTE: most analyses of networks were relegated to htc resources.
 See the README in the `htc` folder for a description of this process. 
These relegated analyses can still be performed locally using the `analysis.jl` script.
Below we describe the analysis scripts that were not used on htc or were ad-hoc analyses.

Here we analyze the properties of the simulated networks, estimated networks, and compare these results. The following files are used here summarize information on the simulated networks:

`sim_sums.R`: This computes some basic information about the network such as the level, number of blobs, class membership (tree-child, tree-based, FU-stable), and confirms the lack of small cycles. Results here are saved to `./data/sim_net_props.csv`

`sim_sums.jl`: This script computes more information about the network, computes specific properties about individual hybridizations, and branch lengths. it will generate the following files and folders:
- - `data/pars_i/hyb_dat.csv`: This contains information about each hybridization found across the networks within the parameterization
- - `data/pars_i/branch_lengths.csv`: This contains the internal branch lengths found across the networks 
- - `data/pars_i/net_j/hybs/H_k/stacked.csv`: This contains the hybrid nodes that are stacked above or below the hybridization
- - `data/pars_i/net_j/hybs/H_k/cycle_nodes.csv`: The nodes that are part of the cycle of the hybridization
- - `data/pars_i/net_j/hybs/H_k/cycle_hybrids.csv`: the nodes that are a part of the cycle of the hybridization and denotes whether they are above or below on the cycle


Note that SNaQ inference depends on the maximum number of hybridizations (hmax) allowed during the search. We ran SNaQ for values of hmax 0-5. This adds another level to our nested dataset, h_no. Generally information is summarizedacross par_no settings (i.e., the csv will be in the par_i folder but information on phy_no, rep_no, and h_no are contained within the csv). it generates the following files in the newly created `summarized_results`.

These files summarize the estimated networks to the true network:

- `net_properties.R`: This computes the tree of blobs for both the estimated and true networks and compares them. It creates the following files:
- - `./summarized_results/par_no/tob.csv` This computes distances (various metrics) between the estimated and true tree of blobs
- - `./summarized_results/par_no/tob/compat.csv` This contains information on whether the individual blobs found within the estimated network are compatible with any of the blobs in the true network
**Note** `net_props_parallel.R` can take advantage of parallelism to perform the same end quicker.

- `analysis.jl`:This script has 3 modes of running depending on the given input argument. **NOTE** that modes **1** and **2** should be redundant as they are computed within the htc resources as a part of the analysis step. However, modes 1 and 2 should still be run as they compile the results into the appropraite csv files. Mode **3**  is needed to compute all neeeded summaries and make figures. This mode was some post-hoc analyses that we wanted to perform after running initial analyses. Alltogether these find the number of 'found' hybrids and clusters in the network, in both compatible and strict senses. It also computes the HWCD between the estimated net and the true network and its subnetworks. It generates three files:
- - **Execution Mode 1 (Main Analysis)**: Iterates across parameter configurations, true networks, and replicates to evaluate topological recovery. It computes hardwired cluster distances (HWCD) against the true network, canonical networks, and fold-stable networks, quantifies hybrid detection performance (true/false positives and negatives), evaluates displayed tree Robinson-Foulds (RF) distances, and integrates goodness-of-fit (GoF) metrics. Generates `clusters.csv` and `found_clades.csv`.
- - **Execution Mode 2 (Squirrel Analysis)**: Transforms networks into semi-network representations to evaluate quartet network consistency scores (`quar_found` and `quar_compat`) using the Squirrel framework. Generates `squirrel.csv`.
- - **Execution Mode 3 (Concordance Factor Analysis)**: Calculates concordance factor (CF) discrepancy metrics—comparing estimated expected CFs to observed CFs, estimated expected CFs to true expected CFs, and true observed CFs—alongside subnetwork topological distances (`subnet_dist`). Generates `cf_dists.csv`.
- - **Outputs of a Run**: Structured summary datasets saved under `../summarized_results/pars_i/` comprising `clusters.csv`, `found_clades.csv`, `squirrel.csv`, and `cf_dists.csv`.

- `regression.R` 
- - Loads network properties from `sim_net_props.csv` and merges them with parameter settings to build a unified analysis of how often hybrid events are found on the estimated networks.
- - **Bootstrap **: Runs 100 bootstrap iterations that downsamples 1 randomly 1 rep for each bootstrap.
- - **Conditional Logistic Regression**: Fits conditional logistic models (`clogit`) to evaluate the impact of network features on exact clade recovery.
- - **Bagged Decision Tree Classification**: Trains bagged decision trees (`treebag`) via 10-fold cross-validation using `caret` to evaluate feature importance and test accuracy.
- - Results  are saved to `../summarized_results/boot_regression.csv` alongside console summaries and `xtable` exports.

`results_summary.R`: filters the data based on a few different criteria and then feeds the filtered datasets to `final-plots.Rmd` to make figures and compute summary statistics.


