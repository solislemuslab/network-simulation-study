library(diversitree)
library(ape)
library(phytools)

set.seed(1)

BiasYuleNumTaxa<-c()
for(i in 1:100){
  yuleTS <- trees(pars=c(1), # in TreeSim or NetSim for Yule: lambda = 1, mu = 0
                type="yule", 
                max.taxa = 25,# in TreeSim or NetSim for Yule: n = 25
                n=100) # in TreeSim or NetSim for Yule: numbsim = 100
  
  # fit a Yule model with diversitree
  yuleTSFits <- lapply(X = yuleTS, 
                       FUN = function(x) diversitree::make.yule(tree = x,sampling.f = 1))
  # fetch the coefficients estimates on the previous lik funs
  yuleTSMLEcoefs <- lapply(X = yuleTSFits, 
                           FUN = function(x) coef(diversitree::find.mle(func = x,x.init = 0.5)))
  #Lambda Bias
  biasMe<-1-mean(unlist(yuleTSMLEcoefs))#Bias
  BiasYuleNumTaxa<-c(BiasYuleNumTaxa,biasMe)
}

#Ploting
hist(BiasYuleNumTaxa)




yuleTS <- trees(c(1), "yule", max.taxa = 25,n=100)
ltt(yuleTS[[1]],show.tree=TRUE)




yuleTS1 <- sim.bd.age(age = 4, numbsim=20, lambda = 1, mu = 0, frac = 1, mrca = TRUE,complete = TRUE, K = 0)

ltt(yuleTS1[[5]],show.tree=TRUE)







