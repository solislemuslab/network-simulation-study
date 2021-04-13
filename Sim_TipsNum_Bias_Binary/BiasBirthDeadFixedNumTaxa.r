

# Lambda bias from birth dead Tree simulation with TreeSim and NetSim with numTaxa=15, numbsim=100
# the simulation was repeated 100 times 

library(NetSim) # simulate bd time-tree networks
library(TreeSim) # simulate bd time-trees
library(ape) # tools for handling phylo objects
library(diversitree) # tools for estimating lambda and mu

#############
BiasLambdaNSn15<-c()

for(i in 1:100){
    numbsim=100
    bdNS <- NetSim::sim.bdh.taxa.ssa(n=15,
                                     numbsim = numbsim,
                                     lambda = 0.9,
                                     mu = 0.5,
                                     nu = 0,
                                     hybprops = c(0.5, 0.5, 0.5),
                                     hyb.inher.fxn = NetSim::make.beta.draw(1, 1),
                                     frac = 0.9,
                                     mrca = TRUE,
                                     complete = FALSE,
                                     stochsampling = FALSE,
                                     hyb.rate.fxn = NULL,
                                     trait.model = NULL)
    ## remove any trees with no taxa
    bdNS <- bdNS[!sapply(X = bdNS, FUN = is.null)]
    ## remove non-phylo elements
    bdNS <- bdNS[sapply(X = bdNS, FUN = is.phylo)]
    #Code to convert not binary Tree in binary
    for(j in 1:length(bdNS)){
        if(is.binary(bdNS[[j]])){
            bdNS[[j]]<-bdNS[[j]]
        } else {
            bdNS[[j]]<-unroot(bdNS[[j]])
        }
    }
    #sapply(X = bdNS, FUN = is.binary)
    # fit a Bd model with diversitree
    bdNSFits <- lapply(X = bdNS, FUN = function(x) diversitree::make.bd(tree = x, sampling.f = 0.6))
    bdNSFits <- lapply(X = bdNS, FUN = function(x) diversitree::make.bd(tree = x, sampling.f = 1))#
    # fetch the coefficients estimates on the previous lik funs
    bdNSMLEcoefs <- lapply(X = bdNSFits, FUN = function(x) coef(diversitree::find.mle(func = x,x.init = c(0.5, 0.5))))
    
    si<-sapply(bdNSMLEcoefs, function(x) unlist(x)["lambda"])
    biansMe<-0.9-mean(si)
    BiasLambdaNSn15<-c(BiasLambdaNSn15,biansMe)
}

BiasLambdaTSn15<-c()

for(i in 1:100){ 
    numbsim = 100
    bdTS <- TreeSim::sim.bd.taxa(n = 15,
                                 numbsim = numbsim,
                                 lambda = 0.9,
                                 mu = 0.5,
                                 frac = 0.9,
                                 complete = FALSE,
                                 stochsampling = TRUE)
    ## remove any trees with no taxa
    bdTS <- bdTS[!sapply(X = bdTS, FUN = is.null)]
    ## remove non-phylo elements
    bdTS <- bdTS[sapply(X = bdTS, FUN = is.phylo)]
    # fit a Bd model with diversitree
    bdTSFits <- lapply(X = bdTS, FUN = function(x) diversitree::make.bd(tree = x, sampling.f = 0.6))
    #bdTSFits <- lapply(X = bdTS, FUN = function(x) diversitree::make.bd(tree = x, sampling.f = 1))#
    # fetch the coefficients estimates on the previous lik funs
    bdTSMLEcoefs <- lapply(X = bdTSFits, FUN = function(x) coef(diversitree::find.mle(func = x,x.init = c(0.5, 0.5))))
    sap<-sapply(bdTSMLEcoefs, function(x) unlist(x)["lambda"])
    #mean(sap)
    #plot(density(sap))
    #abline(v=0.9)
    biansMe<-0.9-mean(sap)
    BiasLambdaTSn15<-c(BiasLambdaTSn15,biansMe)
}

#pdf(file="biasBD.pdf",width=10, height=7)
par(mfrow=c(1,2))
hist(BiasLambdaTSn15,
     main="Bias of Lambda from Birth-death TreeSim\n(100 numbsim repeated 100 times)\n(n=15, lambda=0.9)",
     xlab=expression(paste(lambda, " bias")))
abline(v=0)

hist(BiasLambdaNSn15,
     main="Bias of Lambda from Birth-death NetSim\n(100 numbsim repeated 100 times)\n(n=15, lambda=0.9)",
     xlab=expression(paste(lambda, " bias")))
abline(v=0)
#dev.off()
