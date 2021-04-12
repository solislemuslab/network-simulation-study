
# Lambda bias from yule Tree simulation with TreeSim and NetSim at age = 4 and 6, numbsim=100
# the simulation was repeated 100 times 

library(NetSim) # simulate bd time-tree networks
library(TreeSim) # simulate bd time-trees
library(ape) # tools for handling phylo objects
library(diversitree) # tools for estimating lambda and mu
library(lattice) # stacked plots

numbsim <- 100

# Yule
BiasLamTree4<-c()
for(i in 1:100){
yuleTS <- sim.bd.age(age = 4, numbsim = numbsim, lambda = 1, mu = 0, frac = 1, mrca = TRUE,complete = TRUE, K = 0)
## remove any trees with no taxa
yuleTS <- yuleTS[!sapply(X = yuleTS, FUN = is.null)]
## remove non-phylo elements
yuleTS <- yuleTS[sapply(X = yuleTS, FUN = is.phylo)]
## remove any trees with <10 taxa
yuleTS <- yuleTS[sapply(X = yuleTS, FUN = function(x) length(x$tip.label) >= 10)]
length(yuleTS)
# fit a Yule model with diversitree
yuleTSFits <- lapply(X = yuleTS, FUN = function(x) diversitree::make.yule(tree = x,sampling.f = 1))
# fetch the coefficients estimates on the previous lik funs
yuleTSMLEcoefs <- lapply(X = yuleTSFits, FUN = function(x) coef(diversitree::find.mle(func = x,x.init = 0.5)))
#Lambda Bias
biasYule<-1-mean(unlist(yuleTSMLEcoefs))#Bias
BiasLamTree4<-c(BiasLamTree4,biasYule)
}

BiasLamNetYul4<-c()
for(i in 1:100){
yuleNS <- NetSim::sim.bdh.age(age = 4, numbsim = numbsim, lambda = 1, mu = 0, nu = 0,
hybprops = c(0.5, 0.5, 0.5), hyb.inher.fxn = NetSim::make.beta.draw(1, 1), frac = 1,
mrca = TRUE, complete = TRUE, stochsampling = FALSE, hyb.rate.fxn = NULL, trait.model = NULL)
## remove any trees with no taxa
yuleNS <- yuleNS[!sapply(X = yuleNS, FUN = is.null)]
## remove non-phylo elements
yuleNS <- yuleNS[sapply(X = yuleNS, FUN = is.phylo)]
## remove any trees with <10 taxa
yuleNS <- yuleNS[sapply(X = yuleNS, FUN = function(x) length(x$tip.label) >= 10)]
# fit a Yule model with diversitree
yuleNSFits <- lapply(X = yuleNS, FUN = function(x) diversitree::make.yule(tree = x,
sampling.f = 1))
# fetch the coefficients estimates on the previous lik funs
yuleNSMLEcoefs <- lapply(X = yuleNSFits, FUN = function(x) coef(diversitree::find.mle(func = x,
x.init = 0.5)))
#Lambda Bias
biasYule<-1-mean(unlist(yuleNSMLEcoefs))
BiasLamNetYul4<-c(BiasLamNetYul4,biasYule)
}


# BirthDead

BiasLamTree6<-c()
for(i in 1:100){
yuleTS <- sim.bd.age(age = 6, numbsim = numbsim, lambda = 1, mu = 0, frac = 1, mrca = TRUE,complete = TRUE, K = 0)
## remove any trees with no taxa
yuleTS <- yuleTS[!sapply(X = yuleTS, FUN = is.null)]
## remove non-phylo elements
yuleTS <- yuleTS[sapply(X = yuleTS, FUN = is.phylo)]
## remove any trees with <10 taxa
yuleTS <- yuleTS[sapply(X = yuleTS, FUN = function(x) length(x$tip.label) >= 10)]
length(yuleTS)
# fit a Yule model with diversitree
yuleTSFits <- lapply(X = yuleTS, FUN = function(x) diversitree::make.yule(tree = x,sampling.f = 1))
# fetch the coefficients estimates on the previous lik funs
yuleTSMLEcoefs <- lapply(X = yuleTSFits, FUN = function(x) coef(diversitree::find.mle(func = x,x.init = 0.5)))
# Lambda Bias
biasYule<-1-mean(unlist(yuleTSMLEcoefs))
BiasLamTree6<-c(BiasLamTree6,biasYule)
}




BiasLamNetYul6<-c()
for(i in 1:100){
yuleNS <- NetSim::sim.bdh.age(age = 6, numbsim = numbsim, lambda = 1, mu = 0, nu = 0,
hybprops = c(0.5, 0.5, 0.5), hyb.inher.fxn = NetSim::make.beta.draw(1, 1), frac = 1,
mrca = TRUE, complete = TRUE, stochsampling = FALSE, hyb.rate.fxn = NULL, trait.model = NULL)
## remove any trees with no taxa
yuleNS <- yuleNS[!sapply(X = yuleNS, FUN = is.null)]
## remove non-phylo elements
yuleNS <- yuleNS[sapply(X = yuleNS, FUN = is.phylo)]
## remove any trees with <10 taxa
yuleNS <- yuleNS[sapply(X = yuleNS, FUN = function(x) length(x$tip.label) >= 10)]
# fit a Yule model with diversitree
yuleNSFits <- lapply(X = yuleNS, FUN = function(x) diversitree::make.yule(tree = x,
sampling.f = 1))
# fetch the coefficients estimates on the previous lik funs
yuleNSMLEcoefs <- lapply(X = yuleNSFits, FUN = function(x) coef(diversitree::find.mle(func = x,
x.init = 0.5)))
# Lambda Bias
biasYule<-1-mean(unlist(yuleNSMLEcoefs))
BiasLamNetYul6<-c(BiasLamNetYul6,biasYule)
}


#pdf(file="bias.pdf",width=10, height=7)
par(mfrow=c(2,2))

hist(BiasLamTree4,xlab=expression(paste(lambda, " bias")),main="Yule_TreeSim (100 numbsim repeated 100 times)\n( age=4, lambda=1)")
abline(v=0,col="white",lwd=2)
abline(v=0,lty=3,col="red")

hist(BiasLamNetYul4,xlab=expression(paste(lambda, " bias")),main="Yule_NetSim (100 numbsim repeated 100 times)\n( age=4, lambda=1)")
abline(v=0,col="white",lwd=2)
abline(v=0,lty=3,col="red")

hist(BiasLamTree6,xlab=expression(paste(lambda, " bias")),main="Yule_TreeSim (100 numbsim repeated 100 times)\n( age=6, lambda=1)")
abline(v=0,col="white",lwd=2)
abline(v=0,lty=3,col="red")

hist(BiasLamNetYul6,xlab=expression(paste(lambda, " bias")),main="Yule_NetSim (100 numbsim repeated 100 times)\n( age=6, lambda=1")
abline(v=0,col="white",lwd=2)
abline(v=0,lty=3,col="red")
#dev.off()


