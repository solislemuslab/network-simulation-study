library(diversitree)
library(ape)
library(phytools)

set.seed(10000)

#Yule simulation with trees
# Notice that using the diversitree::trees() function of diversitree, the bias distribution is suitable for Yule tree.

BiasYuleNumTaxa<-c()
for(i in 1:100){
  yuleTS <- diversitree::trees(pars = c(1), # in TreeSim or NetSim for Yule: lambda = 1, mu = 0
                  type = "yule",
				  n = 100,# in TreeSim or NetSim for Yule: numbsim = 100
                  max.taxa = 25,# in TreeSim or NetSim for Yule: n = 25
                  max.t=Inf,#the tree will evolve until ‘max.taxa’ extant taxa are present
				  include.extinct=FALSE#Not include extinct taxa#in TreeSim: Complete=F
				  )
  
  # fit a Yule model with diversitree
  yuleTSFits <- lapply(X = yuleTS,FUN = function(x) diversitree::make.yule(tree = x,sampling.f = 1))
  # fetch the coefficients estimates on the previous lik funs
  yuleTSMLEcoefs <- lapply(X = yuleTSFits,FUN = function(x) coef(diversitree::find.mle(func = x,x.init = 0.5)))
  #Lambda Bias
  biasMe<-1-mean(unlist(yuleTSMLEcoefs))#Bias
  BiasYuleNumTaxa<-c(BiasYuleNumTaxa,biasMe)
}


#Ploting
hist(BiasYuleNumTaxa,
     xlab=expression(paste(lambda, " bias")),
     main="Yule_TreeSim (100 numbsim repeated 100 times)\n( numTaxa=25, lambda=1)")
abline(v=0,col="white",lwd=2)
abline(v=0,lty=3,col="red")




#Birth-Dead simulation, with lambda=1, mu=0.5, numTaxa=15, numbsim = 100

#To estiamte the parameters "optim" method was used in diversitree::find.mle() function (recomend by the author)
#Notice that even using the methodology that diversitree recomend, the bias is negative


BiasLambdaBDNumTaxa<-c()
BiasMuBDNumTaxa<-c()
for(i in 1:100){
  bdTS <- diversitree::trees(pars = c(1,0.5), # in TreeSim or NetSim for BD: lambda = 1, mu = 0.5
                  type = "bd",
				  n = 100,# in TreeSim or NetSim for BD: numbsim = 100
                  max.taxa = 15,# in TreeSim or NetSim for BD: n = 25
                  max.t=Inf,#the tree will evolve until ‘max.taxa’ extant taxa are present
				  include.extinct=F#Not include extinct taxa#in TreeSim: Complete=F
				  )
  
  bdTSFits <- lapply(X = bdTS, FUN = function(x) diversitree::make.bd(tree = x, sampling.f = 1))
  #bdTSFits <- lapply(X = bdTS, FUN = function(x) diversitree::make.bd(tree = x, sampling.f = 1))#
  # fetch the coefficients estimates on the previous lik funs
  bdTSMLEcoefs <- lapply(X = bdTSFits, FUN = function(x) coef(diversitree::find.mle(func = x,x.init = c(0.5, 0.5),method="optim")))
  Lambda <- sapply(bdTSMLEcoefs, function(x) unlist(x)["lambda"])
  Mu <- sapply(bdTSMLEcoefs, function(x) unlist(x)["mu"])
  #Lambda Bias
  biasLambdaMe <- 1-mean(Lambda)#Bias lambda
  biasMuMe <- 0.5-mean(Mu)#Bias mu
  BiasLambdaBDNumTaxa <-c(BiasLambdaBDNumTaxa,biasLambdaMe)
  BiasMuBDNumTaxa <- c(BiasMuBDNumTaxa,biasMuMe)
}

#Ploting
par(mfrow=c(1,2))
hist(BiasLambdaBDNumTaxa,
     xlab=expression(paste(lambda, " bias")),
     main="BD simulation (100 numbsim repeated 100 times)\n( numTaxa=25, lambda=1)")
abline(v=0,col="white",lwd=2)
abline(v=0,lty=3,col="red")

hist(BiasMuBDNumTaxa,
     xlab=expression(paste(mu, " bias")),
     main="Bd Simulation (100 numbsim repeated 100 times)\n( numTaxa=25, lambda=1)")
abline(v=0,col="white",lwd=2)
abline(v=0,lty=3,col="red")

