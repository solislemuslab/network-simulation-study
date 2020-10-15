# Background

- Papers about networks are listed [here](https://github.com/solislemuslab/lab-dynamics/blob/master/lit-review.md)
- Papers that already compared network methods under different simulation scenarios:
    - [Kong2020](https://www.biorxiv.org/content/10.1101/2020.07.27.224022v1). Note that the simulate data on mixtures of trees using ms, not on network directly
    - [Elworth2018](https://arxiv.org/abs/1808.08662). This is a long paper. In section 6.3, they say that they simulate data using ms but on a network directly, not on trees
    - [Hejase2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1277-1)
    - [Kamneva2017](https://journals.sagepub.com/doi/full/10.1177/1176934317691935)
- Missing pieces:
    - a more thorough test of SNaQ when network is not level-1 (SNaQ assumption). I want to simulate networks 
    - a more thorough test on networks that we can expect from real data: depth of hybridization, complexity of hybridizations, combination of short branches + hybridization

- Standard approaches to simulate gene trees under a given network:
    1. ms
    2. [HybridLambda](https://github.com/hybridLambda/hybrid-Lambda) (there is a bug on gene tree branch lengths, but many network methods only use gene tree topologies)

- Standard approaches to simulate a network:
    1. Simulate a tree in R with TreeSim (or ape), and then choose two edges at random, and add a horizontal edge between them. I have code in PhyloNetworks to do this
    2. Using the birth-hybridization model (extension of birth-death to generate trees) in [BEAST2](https://academic.oup.com/mbe/article/35/2/504/4705834)
    3. New upcoming R package [NetSim](https://github.com/jjustison/NetSim) (still unpublished, please don't share)

# Existing scripts
- Useful [scripts](https://github.com/crsl4/InconsistencySpeciesTreeGeneFlow) to simulate gene trees under a network and then estimate the network back with PhyloNet
    - Downsides: scripts are in perl and they use ms which needs a very ugly format for the network. Ideally, we want to write the network in parenthetical format and ms does not accept this format
- [Pipeline](https://github.com/crsl4/PhyloNetworks.jl/wiki) to estimate a network with PhyloNetworks. Good to get started with Julia. Ideally, I would like to do the simulations in Julia to take advantage of all the network functions we already have in PhyloNetworks

# Project plan

At some point, we need to agree on a [project plan](https://github.com/solislemuslab/lab-dynamics/blob/master/project-plan.md). No need to worry about this now.

For now, we focus only on:
- Overall objective: Test the limitations of SNaQ (and other network methods) under a variety of biological scenarios
- Specific aims (need to brainstorm). We want to test performance under:
    - different depths of hybridizations
    - different complexities of hybridizations (level-1 vs not level-1)
    - under short branch lengths
    - ...

# Next steps
(not necessarily in order)

- Read other papers that have done simulations to test networks
- Check existing scripts to see if we can start creating our own simulating scripts for this project
- Think which comparisons we want to make