---
output: html_document
editor_options: 
  chunk_output_type: inline
---
In general we must to have Julia program in tar.gz, and  Julia packages into a project directory in tar.gz. 

Later we have to put this in each job files to run in CHTC with the option: “One job per directory”

## Installation of julia in CHTC

Install Julia and Julia packages using the following instructions:
[https://chtc.cs.wisc.edu/uw-research-computing/julia-jobs.html](https://chtc.cs.wisc.edu/uw-research-computing/julia-jobs.html)


This are the steps for get julia and the project directory as tar.gz:

- Download Julia from web (https://julialang.org/downloads/platform/#linux_and_freebsd) and decompress
  - wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz
  - tar zxvf julia-1.8.5-linux-x86_64.tar.gz

- Install Julia packages in the project directory
  - Create one directory to install Julia packages 
    - mkdir my-project
    - export JULIA_DEPOT_PATH=$PWD/my-project
  - Open Julia inside the directory of Julia packages and install the packages with 
    - julia-1.8.5/bin/julia --project=my-project
      - using Pkg
      - Pkg.add("PhyloNetworks")
      - quit()
  - SAVE INSTALLED PACKAGES FOR LATER JOBS
    - tar -czf pro.tar.gz my-project/
  


## Data organization

HCTC offer many alternatives to parallelize processes in this example we explore the option “One job per directory”, that allow to organize job files into separate directories where each job file directory must to have the same name plus a correlative subindex number that start in zero. In this example we name the directories as chtc0, chtc1, …

In each directory we must to put all the files needed to run the project (programs, libraries and data), and each file have to have the same name in each directory, so to facilitate data management, we put all the data into compress file namely files.tar.gz. To do that we have use **level_2/pipeline/prepare_for_chtc_short.sh** that use the results of **level_2/pipeline/generate_datasets.R** that save the results in **level_2/data**. Finaly to get the data, the script decompress files.tar.gz (toy4.sh). The structure of the file is the next one:

- **files.tar.gz:** Compressed file with all the data needed:
  - **.csv:** File with the concordance factors
  - **.tre:** File with the consensus tree.
  - **.extnewick:** Original network
  - **gt_seed:** Randod seed generated with generate_datasets.R
- **julia-1.6.3-linux-x86_64.tar.gz:** julia porgram.
- **pro.tar.gz:** julia packages project. 
- **network_estimation.jl:** julia script.



## Scripts

- **network_estimation.jl:** Read the concordance factors, the gene tree and the random seed to estimate hybridizations. This script use Distributed package to run in many Threads. First it no run, then Gustavo sugest to use: exeflags=“–project=$(Base.active_project())”: see following link [stackoverflow](https://stackoverflow.com/questions/70792051/package-set-up-not-propagating-to-workers-with-distributed).


  - code sequence.
    - using Distributed
    - addprocs(nthreads, exeflags="--project=$(Base.active_project())")
    - @everywhere using PhyloNetworks
    - @everywhere using Functors


- **toy4.sh:** load and decompress julia (julia-1.6.3-linux-x86_64.tar.gz), julia pakages (pro.tar.gz) and the data (files.tar.gz). Call  the concordance factors, the gene tree files and the random seed and pass it to julia_snaq.jl to estimate hybridizations.

- **toy4.sub:** submit the job.

## Running

### Sh file



```bash
SEED=`cat gt_seed`

julia --project=pro network_estimation.jl *.csv *.tre 0 4 4 `expr $SEED + 1` > snaq_outgroup_h0.outerr 2>&1
julia --project=pro network_estimation.jl *.csv snaq_output_h0.out 1 4 4 `expr $SEED + 3` > snaq_outgroup_h1.outerr 2>&1
julia --project=pro network_estimation.jl *.csv snaq_output_h1.out 2 4 4 `expr $SEED + 5` > snaq_outgroup_h2.outerr 2>&1
julia --project=pro network_estimation.jl *.csv snaq_output_h2.out 3 4 4 `expr $SEED + 7` > snaq_outgroup_h3.outerr 2>&1
```

Line by line:

```bash
julia --project=pro network_estimation.jl *.csv *.tre 0 4 4 `expr $SEED + 1` > snaq_outgroup_h0.outerr 2>&1
```
parameters:

  - `*.csv`: Concordance factor file name.
  - `*.tre`: network file name.
  - 0: subscript of the first output.
  - 4: nthreads.
  - 4: nruns.
  - `expr $SEED + 1`: Random seed.
  
  
  
# prepare_for_chtc_short.sh



```bash
julia --project=pro network_estimation.jl *.csv snaq_output_h0.out 1 4 4 `expr $SEED + 3` > snaq_outgroup_h1.outerr 2>&1
```

parameters:

  - `*.csv`: Concordance factor file name.
  - snaq_output_h0.out: network file name from the first output.
  - 0: subscript of the first output.
  - 4: nthreads.
  - 4: nruns.
  - `expr $SEED + 3`: seed.
