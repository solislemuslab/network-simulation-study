## Installation of julia in CHTC

Install Julia and Julia packages using the following instructions:
[https://chtc.cs.wisc.edu/uw-research-computing/julia-jobs.html](https://chtc.cs.wisc.edu/uw-research-computing/julia-jobs.html)

I did it for julia 1.6, so I will make the same for julia 1.8 and make a short description about that.

In general we must to have julia program in tar.gz, and  julia packages into a project directory in tar.gz. 

Later we have to put this in each job files to run in CHTC with the option: “One job per directory”

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
