
library(NetSim) # simulate bd time-tree networks: https://github.com/jjustison/NetSim
library(geiger)
library(ape) # tools for handling phylo objects

# The goal of this scrip is show how even with lambda = 0.9 and  mu = 0.5 the percentage of extinct trees is much more of 50%



set.seed(300)
numbsim=600
bdNS <- NetSim::sim.bdh.taxa.ssa(n=10, 
                                 numbsim = numbsim, 
                                 lambda = 0.9, 
                                 mu = 0.5, 
                                 nu = 0.03,
                                 hybprops = c(0.5, 0.5, 0.5), 
                                 hyb.inher.fxn = NetSim::make.beta.draw(1, 1), 
                                 frac = 1,
                                 mrca = FALSE, 
                                 complete = TRUE, 
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

# Extract information of reticulate and not reticulate Trees
Ret<-function(x){
NRet<-nrow(x$reticulation)
if(NRet==0){
resp<-"Not Reticulation"
}else if(NRet>0){
resp<-"Reticulation Tree"
}
return(resp)
}

RetChar<-sapply(X = bdNS2, FUN = Ret)

# Percentaje of tree types
AllTreeType<-c(ExtTip1Char,RetChar)
PerTree<-(table(AllTreeType)/length(AllTreeType))*100

#Making plot
barplot(PerTree,main="Percentage of tree types",ylim=c(0,60),col=c("#ff6f69","#88d8b0","#88d8b0"))

