library(NetSim) # simulate bd time-tree networks: https://github.com/jjustison/NetSim
library(TreeSim) # simulate bd time-trees
library(ape) # tools for handling phylo objects
library(diversitree) # tools for estimating lambda and mu
library(lattice) # stacked plots
set.seed(1000)
numbsim <- 100

# Check in which situations the `.age` functions in `TreeSim` and `NetSim` produce binary trees. The whole thing should be controlled through the argument `mrca` in both packages. Also, it is possible that malformed trees come from short-lived processes where the whole tree becomes extinct or just fails to speciate; this can be explored using three arbitrary levels of time: 4, 8, and 10 time units.

# In addition, the function `ape::is.binary` tests whether the tree is binary so the whole trick is to test with mrca = TRUE and FALSE for each packages (i.e., 12 simulations), then apply the function `ape::is.binary` to each tree in a given collection and summarize the number of binary trees (if any).

# PLEASE NOTE THAT THIS CODE TAKES ABOUT AN HOUR OR SO FOR RUNNING. AN .RDATA FILE IS AVAILABLE UPON REQUEST AS .GITIGNORE WON'T ALLOW TO ADD SUCH FILES.

#repository of simulation results: package-mrca-age
mrca <- c(TRUE, FALSE)
age <- c(4, 6, 8)
simNames <- apply(X = expand.grid(mrca,
                                  age),
                  MARGIN = 1,
                  FUN = function(x) paste(x, collapse = "-"))
simTreeSim <- list()
length(simTreeSim) <- length(simNames)
#names(simTreeSim) <- simNames
simNetSim <- list()
length(simNetSim) <- length(simNames)
#names(simNetSim) <- simNames

# Begin simulating trees for `TreeSim`:

# using TreeSim
nIter <- 1
for (i in mrca) {
    for (j in age) {
        cat("Creating sim with TreeSim, ", "mrca = ", i, " age = ", j, "\n", sep = "")
        simTreeSim[nIter] <- list(eval(call("sim.bd.age",
                                       age = j,
                                       numbsim = numbsim,
                                       lambda = 1,
                                       mu = 0,
                                       frac = 1,
                                       mrca = FALSE,
                                       complete = TRUE,
                                       K = 0)))
        simNames[nIter] <- paste("TS", i, j, sep = "_")
        print(simNames[nIter])
        nIter <- nIter + 1
    }
}
names(simTreeSim) <- simNames

# Now simulating trees for `NetSim`:

# using NetSim
nIter <- 1
for (i in mrca) {
    for (j in age) {
        cat("Creating sim with NetSim, ", "mrca = ", i, " age = ", j, "\n", sep = "")
        simNetSim[[nIter]] <- eval(call("sim.bdh.age",
                                      age = j,
                                      numbsim = numbsim,
                                      lambda = 1,
                                      mu = 0,
                                      nu = 0,
                                      hybprops = c(0.5, 0.5, 0.5),
                                      hyb.inher.fxn = NetSim::make.beta.draw(1,1),
                                      frac = 1,
                                      mrca = i,
                                      complete = TRUE,
                                      stochsampling = FALSE,
                                      hyb.rate.fxn = NULL,
                                      trait.model = NULL))
        simNames[nIter] <- paste("NS", i, j, sep = "_")
        print(simNames[nIter])
        nIter <- nIter + 1
    }
}
names(simNetSim) <- simNames

# save precious time saving the simulation results to a file until my system comes back
save(simTreeSim, simNetSim, file = "TreeSim_NetSim.Rdata")
