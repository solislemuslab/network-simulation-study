##Libraries
library(NetSim) # simulate bd time-tree networks: https://github.com/jjustison/NetSim
library(geiger)
library(ape) # tools for handling phylo objects


#The goald of this scrip is show that even using mrca=TRUE in bd.age() and bdh.taxa.ssa() it produce binary and not binary trees, and as we expected mrca=FALSE produce all not binary trees



# Birth-Dead Tree sim.bdh.age

## sim.bdh.age (mrca=TRUE)
set.seed(20)
numbsim=300

bdNS <- NetSim::sim.bdh.age(age = 6,
                            numbsim = numbsim,
                            lambda = 0.9,
                            mu = 0.5,
                            nu = 0,
                            hybprops = c(0.5, 0.5, 0.5),
                            hyb.inher.fxn = NetSim::make.beta.draw(1, 1),
                            frac = 1,
                            mrca = TRUE,
                            complete = FALSE,
                            stochsampling = FALSE,
                            hyb.rate.fxn = NULL,
                            trait.model = NULL)

## remove any trees with no taxa
bdNS <- bdNS[!sapply(X = bdNS, FUN = is.null)]

# Extract information of extinct Trees and Trees with one tip
ExtTip1<-unlist(bdNS[which(!sapply(X = bdNS, FUN = is.phylo))])
ExtTip1Char<-as.character(ExtTip1)
ExtTip1Char[ExtTip1Char=="0"]<-"Extinct"
ExtTip1Char[ExtTip1Char=="1"]<-"OneTip"

#remove extinct Trees and Trees with one tip
bdNS2 <- bdNS[sapply(X = bdNS, FUN = is.phylo)]#birth dead tree with NetSim with value 0 or 1

# Extract information of Binary and not binary Trees
Bin<-sapply(X = bdNS2, FUN = is.binary)
BinCh<-as.character(Bin)
BinCh[BinCh=="TRUE"]<-"Binary"
BinCh[BinCh=="FALSE"]<-"NotBinary"


# Percentaje of tree types
AllTreeType<-c(ExtTip1Char,BinCh)
AllTreeType<-factor(AllTreeType,levels=c("Extinct","OneTip","Binary","NotBinary"))
PerTree<-(table(AllTreeType)/length(AllTreeType))*100

#Making plot
#As you can see, even with mrca=TRUE, the simulation produce binary and not binary trees, moreover not Binary Tree is more than 50% of the simulated trees.
par(mfrow=c(1,3),oma=c(0,0,5,0))
barplot(PerTree,main="Percentage of tree types",ylim=c(0,50))
plot(bdNS2[[which(Bin)[1]]],main="First Binary Tree")
plot(bdNS2[[which(!Bin)[1]]],main="First Not Binary Tree")
title("300 Birth-Dead Tree \n (age=6,lambda=0.9,mu=0.5,nu=0,hybprops = c(0.5, 0.5, 0.5))\n (mrca=TRUE)", outer=TRUE,cex.main=1.5) 




## sim.bdh.age (mrca=FALSE)
set.seed(40)
numbsim=300

bdNS <- NetSim::sim.bdh.age(age = 6,
                            numbsim = numbsim,
                            lambda = 0.9,
                            mu = 0.5,
                            nu = 0,
                            hybprops = c(0.5, 0.5, 0.5),
                            hyb.inher.fxn = NetSim::make.beta.draw(1, 1),
                            frac = 1,
                            mrca = TRUE,
                            complete = FALSE,
                            stochsampling = FALSE,
                            hyb.rate.fxn = NULL,
                            trait.model = NULL)

## remove any trees with no taxa
bdNS <- bdNS[!sapply(X = bdNS, FUN = is.null)]

# Extract information of extinct Trees and Trees with one tip
ExtTip1<-unlist(bdNS[which(!sapply(X = bdNS, FUN = is.phylo))])
ExtTip1Char<-as.character(ExtTip1)
ExtTip1Char[ExtTip1Char=="0"]<-"Extinct"
ExtTip1Char[ExtTip1Char=="1"]<-"OneTip"

#remove extinct Trees and Trees with one tip
bdNS2 <- bdNS[sapply(X = bdNS, FUN = is.phylo)]#birth dead tree with NetSim with value 0 or 1

# Extract information of Binary and not binary Trees
Bin<-sapply(X = bdNS2, FUN = is.binary)
BinCh<-as.character(Bin)
BinCh[BinCh=="TRUE"]<-"Binary"
BinCh[BinCh=="FALSE"]<-"NotBinary"


# Percentaje of tree types
AllTreeType<-c(ExtTip1Char,BinCh)
AllTreeType<-factor(AllTreeType,levels=c("Extinct","OneTip","Binary","NotBinary"))
PerTree<-(table(AllTreeType)/length(AllTreeType))*100

#Making plot
par(mfrow=c(1,2),oma=c(0,0,5,0))
barplot(PerTree,main="Percentage of tree types",ylim=c(0,50))
plot(bdNS2[[which(!Bin)[1]]],main="First Not Binary Tree")
title("300 Birth-Dead Tree \n (age=6,lambda=0.9,mu=0.5,nu=0,hybprops = c(0.5, 0.5, 0.5))\n (mrca=TRUE)", outer=TRUE,cex.main=1.5) 





# Birth-Dead Tree sim.bdh.taxa.ssa

## sim.bdh.taxa.ssa (mrca=TRUE)
set.seed(20)
numbsim=300

bdNS <- NetSim::sim.bdh.taxa.ssa(n = 10,
                            numbsim = numbsim,
                            lambda = 0.9,
                            mu = 0.5,
                            nu = 0,
                            hybprops = c(0.5, 0.5, 0.5),
                            hyb.inher.fxn = NetSim::make.beta.draw(1, 1),
                            frac = 1,
                            mrca = TRUE,
                            complete = FALSE,
                            stochsampling = FALSE,
                            hyb.rate.fxn = NULL,
                            trait.model = NULL)


## remove any trees with no taxa
bdNS <- bdNS[!sapply(X = bdNS, FUN = is.null)]

# Extract information of extinct Trees and Trees with one tip
ExtTip1<-unlist(bdNS[which(!sapply(X = bdNS, FUN = is.phylo))])
ExtTip1Char<-as.character(ExtTip1)
ExtTip1Char[ExtTip1Char=="0"]<-"Extinct"
ExtTip1Char[ExtTip1Char=="1"]<-"OneTip"

#remove extinct Trees and Trees with one tip
bdNS2 <- bdNS[sapply(X = bdNS, FUN = is.phylo)]#birth dead tree with NetSim with value 0 or 1

# Extract information of Binary and not binary Trees
Bin<-sapply(X = bdNS2, FUN = is.binary)
BinCh<-as.character(Bin)
BinCh[BinCh=="TRUE"]<-"Binary"
BinCh[BinCh=="FALSE"]<-"NotBinary"


# Percentaje of tree types
AllTreeType<-c(ExtTip1Char,BinCh)
AllTreeType<-factor(AllTreeType,levels=c("Extinct","OneTip","Binary","NotBinary"))
PerTree<-(table(AllTreeType)/length(AllTreeType))*100

#Making plot
par(mfrow=c(1,3),oma=c(0,0,5,0))
barplot(PerTree,main="Percentage of tree types",ylim=c(0,50))
plot(bdNS2[[which(Bin)[1]]],main="First Binary Tree")
plot(bdNS2[[which(!Bin)[1]]],main="First Not Binary Tree")
title("300 Birth-Dead Tree \n (n=10,lambda=0.9,mu=0.5,nu=0,hybprops = c(0.5, 0.5, 0.5))\n (mrca=TRUE)", outer=TRUE,cex.main=1.5) 




## sim.bdh.taxa.ssa (mrca=FALSE)
set.seed(40)
numbsim=300

bdNS <- NetSim::sim.bdh.taxa.ssa(n = 10,
                            numbsim = numbsim,
                            lambda = 0.9,
                            mu = 0.5,
                            nu = 0,
                            hybprops = c(0.5, 0.5, 0.5),
                            hyb.inher.fxn = NetSim::make.beta.draw(1, 1),
                            frac = 1,
                            mrca = TRUE,
                            complete = FALSE,
                            stochsampling = FALSE,
                            hyb.rate.fxn = NULL,
                            trait.model = NULL)

## remove any trees with no taxa
bdNS <- bdNS[!sapply(X = bdNS, FUN = is.null)]

# Extract information of extinct Trees and Trees with one tip
ExtTip1<-unlist(bdNS[which(!sapply(X = bdNS, FUN = is.phylo))])
ExtTip1Char<-as.character(ExtTip1)
ExtTip1Char[ExtTip1Char=="0"]<-"Extinct"
ExtTip1Char[ExtTip1Char=="1"]<-"OneTip"

#remove extinct Trees and Trees with one tip
bdNS2 <- bdNS[sapply(X = bdNS, FUN = is.phylo)]#birth dead tree with NetSim with value 0 or 1

# Extract information of Binary and not binary Trees
Bin<-sapply(X = bdNS2, FUN = is.binary)
BinCh<-as.character(Bin)
BinCh[BinCh=="TRUE"]<-"Binary"
BinCh[BinCh=="FALSE"]<-"NotBinary"


# Percentaje of tree types
AllTreeType<-c(ExtTip1Char,BinCh)
AllTreeType<-factor(AllTreeType,levels=c("Extinct","OneTip","Binary","NotBinary"))
PerTree<-(table(AllTreeType)/length(AllTreeType))*100

#Making plot
par(mfrow=c(1,2),oma=c(0,0,5,0))
barplot(PerTree,main="Percentage of tree types",ylim=c(0,50))
plot(bdNS2[[which(!Bin)[1]]],main="First Not Binary Tree")
title("300 Birth-Dead Tree \n (n=10,lambda=0.9,mu=0.5,nu=0,hybprops = c(0.5, 0.5, 0.5))\n (mrca=TRUE)", outer=TRUE,cex.main=1.5) 



