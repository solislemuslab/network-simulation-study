---
editor_options: 
  markdown: 
    wrap: 72
---

# Pipeline for assessing the behavior of snaq on level-2 networks

The tasks accomplished by each of the scripts is described below..
Please note that execution should be carried out as indicated.

## `Species_Gene_Networks.sh`

### What the script does.

-   Create all the path structures to store all the results in an
    orderly way.

-   Simulate species Networks with SimSpeciesNetworks.r code.

-   Select only the species Networks with no single reticulation (In the
    future we are going to consider all the species Networks).

-   Save the spesie Networks in the temporary path "Output/RNetworks".

-   Read the simulated species Networks with julia, save it as julia
    format in the temporary path "Output/RJuliaNetworks" and in the
    Hybrid lambda format in the temporary path "Output/JuliaForHybLam".

-   Copy and paste the species Networks in Hybrid lambda format to "src"
    of hybrid lambda path, run Hybrid lambda and move all the results
    into the temporary path "Output/TemporaryHyLambdaStorage".

-   Evaluate the Hybrid lambda results to know if the specie networks
    are considered for Hybrid lambda as Ultrametric or not and save the
    previous results in the path "Output/UltraMetric" and
    "Output/NonUltraMetric" as appropriate using the script
    UltrametricEvaluation.r.

-   As we use hybrid lambda, at the moment we will only work with the
    ultrametric species networks.

-   Generate "-n" species tree from each ultrametric species networks
    (storage in "Output/UltraMetric/Ult_JuliaForHybLam") and save it in
    the path "Output/UltraMetric/GeneTrees" and also copy in the path
    "Output/SnaqOut/GeneTrees" to use later with snaq to reconstruct the
    reticulations.

**What are its inputs.**

we have to specify the hybrid lambda src path in the bash script (line
20) Then we have to write the follow parameters:

-   t: It manages two R parameters separated with ";" and without space.

-   numbsim1: Number of Net works.

-   n1: Number of tips.

-   s: It manages four R parameters separated with ";" and without
    space.

-   lambda:

    -   mu=0: extinction parameter, we use 0 (Yule model).

    -   nu: parameter of hybridization.

    -   hybprops: hybridation rates.

-   n: Number of gene trees (hybrid lambda parameter).

-   d: seed of hybrid lambda to generate gene trees.

Example: bash scriptPru.sh -t "numbsim1=500;n1=8" -s
"lambda=0.9;mu=0;nu=0.09;hybprops=c(0.5,0.25,0.25)" -n 500 -d 2

**What its outputs.**

Species networks with no single reticulation. Gene trees from the
ultrametric species networks determined by hybrid lambda.

Note: The number of ultrametric species networks in this moment depends
on the R parameters (number of species networks with no single
reticulation) and the ultrametric hybrid lambda criterion. We can use
the "while" function and repeat the process until we get a specific
number of ultrametric species networks.

**Which dependencies is has.**

-   R libraries:

    -   SiPhyNetwork

    -   geiger

-   julia libraries:

    -   PhyloNetworks

    -   PhyloPlots

    -   RCall

    -   CSV

-   hybrid-Lambda-0.6.2-beta

**Which scripts it runs automatically.**

The following r libraries: - SimSpeciesNetworks.r -
Net_HybridLamb_Format.jl - UltrametricEvaluation.r

## `SimSpeciesNetworks.r`

**What the script does.**

Simulates specific number of networks under the following parameters:

-   numbsim1: Number of Net works  
-   n1: Number of tips  
-   lambda:  
-   mu=0: extinction parameter, we use 0 (Yule model)  
-   nu: parameter of hybridization  
-   hybprops: hybridation rates. Then it evaluates the networks and
    selects the networks with at least one no single reticulation.

**What are its inputs.**

The R parameters

**What its outputs.**

Species networks with at least one no single reticulation.

**Which dependencies is has.**

-   SiPhyNetwork
-   geiger

**Which scripts it runs automatically.** Nothing

## `Net_HybridLamb_Format.jl`

**What the script does.**

Reads the Species networks, saves the parenthetical format as julia and
hybrid lambda format in the temporary path "Output/RJuliaNetworks" and
"Output/JuliaForHybLam" respectively.

**What are its inputs.**

Species networks with at least one no single reticulation.

**What its outputs.**

Parenthetical format as julia and hybrid lambda format.

**Which dependencies is has.**

Julia libraries:

-   PhyloNetworks
-   PhyloPlots
-   RCall

**Which scripts it runs automatically.**

Nothing

## `UltrametricEvaluation.r`

**What the script does.**

Reads the hybrid lambda output (outerrHybUlt_num.txt) where say if is or
not Ultrametric and then puts all the files from the temporary path in
the ultrametric or no ultrametric path as appropriate and finally
deletes the temporary path.

**What are its inputs.**

hybrid lambda output (outerrHybUlt_num.txt)

**What its outputs.**

Nothing, only organize the files in two main groups (Ultrametric and not
Ultrametric)

**Which dependencies is has.**

Nothing

**Which scripts it runs automatically.**

Nothing

## `SnaqHybridzation.jl`

**What the script does.**

-   Read the gene trees from the path "Output/SnaqOut/GeneTrees", get a
    consensus tree with Quartet MaxCut (namTree) and the concordance
    factors (raxmlCF), finally save that in the path
    "Output/SnaqOut/CFactors". Notice that this process is independent,
    so we can parallelize that.
-   Runs snaq to find reticulations using the concordance factors and
    the consensus tree, finally save the Network in the path
    "Output/SnaqOut/SnaqNet".

**What are its inputs.**

Gene trees

**What its outputs.**

-   Consensus tree from the gene trees
-   Concordance factors from the gene trees

**Which dependencies is has.**

At the moment the path, we can improve that, but the path of Quartet
MaxCut is needed.

**Which scripts it runs automatically.**

Nothing
