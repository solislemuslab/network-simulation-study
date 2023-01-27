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
    - tar -czf my-project.tar.gz my-project/
  


## Data organization

HCTC offer many alternatives to parallelize processes in this example we explore the option “One job per directory”, that allow to organize job files into separate directories where each job file directory must to have the same name plus a correlative subindex number that start in zero. In this example we name the directories as julia_snaq0, julia_snaq1, …

In each directories we must to put the inputs whit the same name, in this example we have the following inputs per job:

- **tableCF.csv:** for the concordance factors.
- **raxmltrees.tre.txt:** for the gene tree and consensus tree.
- **julia-1.6.3-linux-x86_64.tar.gz:** julia porgram.
- **my-julia-project.tar.gz:** julia packages.
- **julia_snaq.jl:** julia script.



## Scripts

- **julia_snaq.jl:** read the concordance factors and the gene tree to estimate hybridizations 

- **julia_snaq.sh:** load julia and julia pakages stored in tar.gz. Call  the concordance factors and the gene tree files and pass it to julia_snaq.jl to estimate hybridizations.

- **julia_snaq.sub:** submit the job.
