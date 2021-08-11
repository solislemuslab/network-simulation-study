---
title: "Level-2 networks sometimes cause segfault or warnings about ultrametricity"
author: Gustavo A. Ballen, Carlos Acosta, and Claudia Solís-Lemus
date: 14-jul-2021
---

Hello,

We have been simulating phylogenetic networks with `SiPhyNetwork` in R, and are trying to use `hybrid-Lambda` to simulate gene trees on these. However, when using some of our level-2 networks as input for hybrid-Lambda we get (1) segmentation faults, (2) a warning on lack of ultrametricity, or (3) both at the same time.

We have used the approach below for network generation:

```{R}
### script.R

library(SiPhyNetwork)
library(ape)
set.seed(10)

# simulate networks under the birth-death hybridization process
numbsim = 500
n1 <- 15
bdNS <- SiPhyNetwork::sim.bdh.taxa.ssa(
    n = n1,
    numbsim = numbsim,
    lambda = 0.9,
    mu = 0,
    nu = 0.03,
    hybprops = c(0.5, 0.25, 0.25),
    hyb.inher.fxn = SiPhyNetwork::make.beta.draw(1, 1),
    frac = 1,
    mrca = FALSE,
    complete = TRUE,
    stochsampling = FALSE,
    hyb.rate.fxn = NULL,
    trait.model = NULL
    )

# here some code testing whether the networks are level-1 or not
out <- # indices to networks that are not level-1

# check that all non-level-1 networks are ultrametric using ape
if (sum(sapply(X = bdNS[out], FUN = ape::is.ultrametric) == length(out))) {
    cat("All networks to be written to nets_not_level1.extnex are ultrametric\n")
}

# write all the networks above level-1 to a file (all of them are)
SiPhyNetwork::write.net(net = bdNS[out], file = "nets_not_level1.extnex")
```

The content of `net_not_level1.extnex` is found attached to this issue.

Then, we use `PhyloNetworks` for conversion of format between extended newick and hybridlambda and write them to the file `nets_not_level1.hybridlambda` (also available in the attachment):

```{julia}
### script.jl

using PhyloNetworks
nets_not_level1 = readMultiTopology("nets_not_level1.extnex")

# bulk format conversion of networks read in nets_not_level1
open("nets_not_level1.hybridlambda", "w") do io    
    for i = 1:length(nets_not_level1)
        hyblambFormatted = hybridlambdaformat(nets_not_level1[i])
        write(io, "$hyblambFormatted\n")
        flush(io)
    end
end
```

And finally iterate over each network to simulate a set of three loci with some arbitrary parameters:

```{bash}
### script.sh
#!/usr/bin/bash

# run the r script
Rscript script.r

# run the julia script
julia script.jl

# iterate over the content of nets_not_level1.hybridlambda and
# run hybrid-Lambda on each line in order to check issues
# with ultrametricity
while read line
do
    let LNUM+=1
    echo "Processing network in line $LNUM:"
    hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
done < nets_not_level1.hybridlambda
```

The above is run redirecting std out and err to a file (also attached):

```{bash}
./script.sh > outerr.txt 2>&1
```

The whole procedure returns 63 networks stored in `nets_not_level1.hybridlambda`, some of which show the problems already mentioned:

- segfault networks: 5, 6, 14, 16, 20, 27, 28, 29, 32, 36, 40, 43, 44, 50,  56, 59, 63
- ultrametric and segfault networks: 34, 52, 57, 58, 61

We would like to know if we are missing something on the way hybrid-Lambda should be run, e.g., whether level-2 networks are not guaranteed to work.

Also, we would like to understand better the warning about lack of ultrametricity, and whether it is safe to ignore it due to the fact that the simulated network is by definition ultrametric, and tested positive with `ape`.

Finally, we would like to understand why the segfault is happening sometimes, and whether it can be solved somehow.

Thanks in advance for your time and feedback.

Best wishes,

Gustavo, Carlos, and Claudia.

# reply to Joe's message

Hello Joe,

Thanks for taking a look into this. I checked out into release `v0.6.2-beta` (whose version returns v0.6.1-beta after building, though):

```bash
hybrid-Lambda/src$ git checkout v0.6.2-beta
HEAD is now at de99ef3 update
hybrid-Lambda/src$ git status
HEAD detached at v0.6.2-beta
hybrid-Lambda/src$ make
hybrid-Lambda/src$ ./hybrid-Lambda
hybrid-Lambda v0.6.1-beta
Usage:
...
```
I ran again the script and all previous instances of segfaults are now gone. However, some of these are instead reported as not ultrametric now, or a couple instances which were initially reported as ultrametric are no longer so:
```diff
diff --git a/outerr.txt b/outerr.txt
index e076b82..5e5f413 100644
--- a/BashScript_Network_NotLevelOne/outerr.txt
+++ b/BashScript_Network_NotLevelOne/outerr.txt
@@ -42,11 +42,19 @@ level2_2_coal_unit
 Processing network in line 5:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
-./script.sh: line 23:  8781 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 6:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
-./script.sh: line 23:  8784 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+WARNING! NOT ULTRAMETRIC!!!
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 7:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -90,7 +98,11 @@ Random seed: 2
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 14:
-./script.sh: line 23:  8802 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 15:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -98,11 +110,19 @@ Random seed: 2
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 16:
-./script.sh: line 23:  8812 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 17:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
+WARNING! NOT ULTRAMETRIC!!!
 Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 18:
@@ -118,7 +138,15 @@ Random seed: 2
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 20:
-./script.sh: line 23:  8824 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+WARNING! NOT ULTRAMETRIC!!!
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 21:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -140,7 +168,11 @@ level2_2_coal_unit
 Processing network in line 24:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
+WARNING! NOT ULTRAMETRIC!!!
 Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 25:
@@ -156,13 +188,31 @@ Random seed: 2
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 27:
-./script.sh: line 23:  8838 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
4:35
Processing network in line 28:
-./script.sh: line 23:  8844 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+WARNING! NOT ULTRAMETRIC!!!
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 29:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
-./script.sh: line 23:  8852 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+WARNING! NOT ULTRAMETRIC!!!
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 30:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -176,7 +226,11 @@ Random seed: 2
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 32:
-./script.sh: line 23:  8865 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 33:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -187,7 +241,12 @@ Processing network in line 34:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
 WARNING! NOT ULTRAMETRIC!!!
-./script.sh: line 23:  8874 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 35:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -197,7 +256,9 @@ level2_2_coal_unit
 Processing network in line 36:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
-./script.sh: line 23:  8883 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 37:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -219,7 +280,9 @@ level2_2_coal_unit
 Processing network in line 40:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
-./script.sh: line 23:  8894 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 41:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -235,9 +298,15 @@ level2_2_coal_unit
 Processing network in line 43:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
-./script.sh: line 23:  8904 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 44:
-./script.sh: line 23:  8909 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 45:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -247,7 +316,11 @@ level2_2_coal_unit
 Processing network in line 46:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
+WARNING! NOT ULTRAMETRIC!!!
 Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 47:
@@ -269,7 +342,15 @@ Random seed: 2
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 50:
-./script.sh: line 23:  8925 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+WARNING! NOT ULTRAMETRIC!!!
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 51:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -279,8 +360,9 @@ level2_2_coal_unit
 Processing network in line 52:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
-WARNING! NOT ULTRAMETRIC!!!
-./script.sh: line 23:  8934 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 53:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -290,7 +372,11 @@ level2_2_coal_unit
 Processing network in line 54:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
+WARNING! NOT ULTRAMETRIC!!!
 Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 55:
@@ -302,19 +388,35 @@ level2_2_coal_unit
 Processing network in line 56:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
-./script.sh: line 23:  8947 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 57:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
 WARNING! NOT ULTRAMETRIC!!!
-./script.sh: line 23:  8950 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 58:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
 WARNING! NOT ULTRAMETRIC!!!
-./script.sh: line 23:  8958 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 59:
-./script.sh: line 23:  8966 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 60:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
@@ -325,15 +427,28 @@ Processing network in line 61:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
 WARNING! NOT ULTRAMETRIC!!!
-./script.sh: line 23:  8978 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 62:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
+WARNING! NOT ULTRAMETRIC!!!
 Random seed: 2
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
+WARNING: Gene tree is not ultrametric
 Produced gene tree files:
 level2_2_coal_unit
 Processing network in line 63:
-./script.sh: line 23:  8988 Segmentation fault      (core dumped) hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
+Default Kingman coalescent on all branches.
+Default population size of 10000 on all branches.
+Random seed: 2
+Produced gene tree files:
+level2_2_coal_unit
 Processing network in line 64:
 Default Kingman coalescent on all branches.
 Default population size of 10000 on all branches.
```
Hope these results help to solve the bug. We are still trying to understand what triggers the warning, and whether it is safe to use hybrid-Lambda in the presence of such lack of ultrametricity. We may discard these not-ultrametric networks if it is better to do so, but given that these were simulated and in principle should all be ultrametric, we are still puzzled by this warning.
By looking at the code of `net.cpp` ([L323](https://github.com/hybridLambda/hybrid-Lambda/blob/e782a53947ba0aab3967edcc2e53bb9466a3304e/src/net.cpp#L323)) we find that the function `check_isUltrametric` uses `0.000001` in square difference of two paths as a threshold for considering them different, and hence the network as not ultrametric. That makes absolute differences in path length greater than 0.001 to fail the test. Is there any reason for choosing such particular threshold and not perhaps a bigger value? That could help us understand whether the simulation of the networks might produce misleading results.
Thanks again,
Gustavo
