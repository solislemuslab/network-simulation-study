#Tips_Number distribution of yule tree with NetSim and TreeSim with numbsim=500#

#Diversitree need number of tips >=10 to estimate parameteres
#so we have to know the percentage of data remaining
#Fortunately with age=4 we keep almost 100% of the data when we use mrca=TRUE
#but when we work with mrca=FALSE we kep almos 90% of the data.
# In adition to we can see that the distribution is asimetric


library(NetSim) # simulate bd time-tree networks
library(TreeSim) # simulate bd time-trees
library(ape) # tools for handling phylo objects
library(diversitree) # tools for estimating lambda and mu

numbsim <- 500

#TreeSim age=4
yuleTS <- TreeSim::sim.bd.age(age = 4,
                              numbsim = numbsim,
                              lambda = 1,
                              mu = 0,
                              frac = 1,
                              mrca = FALSE,
                              complete = TRUE,
                              K = 0)
## remove any trees with no taxa
yuleTS <- yuleTS[!sapply(X = yuleTS, FUN = is.null)]
## remove non-phylo elements
yuleTS <- yuleTS[sapply(X = yuleTS, FUN = is.phylo)]
## Tip Number 
TipYulTS4 <- unlist(sapply(X = yuleTS, FUN = function(x) length(x$tip.label)))
## Percentage of tree with tip number >=10
PorTip10YuTS4<-round(mean(TipYulTS4>=10)*100,3)

#NetSim age=4
yuleNS <- NetSim::sim.bdh.age(age = 4,
                              numbsim = numbsim,
                              lambda = 1,
                              mu = 0,
                              nu = 0,
                              hybprops = c(0.5, 0.5, 0.5),
                              hyb.inher.fxn = NetSim::make.beta.draw(1, 1),
                              frac = 1,
                              mrca = TRUE,
                              complete = TRUE,
                              stochsampling = FALSE,
                              hyb.rate.fxn = NULL,
                              trait.model = NULL)
## remove any trees with no taxa
yuleNS <- yuleNS[!sapply(X = yuleNS, FUN = is.null)]
## remove non-phylo elements
yuleNS <- yuleNS[sapply(X = yuleNS, FUN = is.phylo)]
## Tip Number 
TipYulNS4  <- unlist(sapply(X = yuleNS, FUN = function(x) length(x$tip.label)))
## Percentage of tree with tip number >=10
PorTip10YuNS4 <- round(mean(TipYulNS4>=10)*100,3)

#TreeSim age=6
yuleTS <- TreeSim::sim.bd.age(age = 6,
                              numbsim = numbsim,
                              lambda = 1,
                              mu = 0,
                              frac = 1,
                              mrca = TRUE,
                              complete = TRUE,
                              K = 0)
## remove any trees with no taxa
yuleTS <- yuleTS[!sapply(X = yuleTS, FUN = is.null)]
## remove non-phylo elements
yuleTS <- yuleTS[sapply(X = yuleTS, FUN = is.phylo)]
## Tip Number 
TipYulTS7 <- unlist(sapply(X = yuleTS, FUN = function(x) length(x$tip.label)))
## Percentage of tree with tip number >=10
PorTip10YuTS7<-round(mean(TipYulTS7>=10)*100,3)

#NetSim age=6
yuleNS <- NetSim::sim.bdh.age(age = 6,
                              numbsim = numbsim,
                              lambda = 1,
                              mu = 0,
                              nu = 0,
                              hybprops = c(0.5, 0.5, 0.5),
                              hyb.inher.fxn = NetSim::make.beta.draw(1, 1),
                              frac = 1,
                              mrca = TRUE,
                              complete = TRUE,
                              stochsampling = FALSE,
                              hyb.rate.fxn = NULL,
                              trait.model = NULL)
## remove any trees with no taxa
yuleNS <- yuleNS[!sapply(X = yuleNS, FUN = is.null)]
## remove non-phylo elements
yuleNS <- yuleNS[sapply(X = yuleNS, FUN = is.phylo)]
## Tip Number 
TipYulNS7  <- unlist(sapply(X = yuleNS, FUN = function(x) length(x$tip.label)))
## Percentage of tree with tip number >=10
PorTip10YuNS7 <- round(mean(TipYulNS7>=10)*100,3)

#pdf(file="Ki.pdf",width=10, height=7)
par(mfrow=c(2,2))
h <- hist(TipYulTS4,
          main="Tip Number from Yule simulations TreeSim\n( numSimb=500, age=4, lambda=1)",
          xlab="Tips")
text(x=quantile(h$breaks,0.7),
     y=quantile(h$counts,0.7),
     labels=paste("Tips>=10: ",PorTip10YuTS4," %"),
     cex=1.3)

h <- hist(TipYulNS4,
          main="Tip Number from Yule simulations NetSim\n( numSimb=500, age=4, lambda=1)",
          xlab="Tips")
text(x=quantile(h$breaks,0.7),
     y=quantile(h$counts,0.7),
     labels=paste("Tips>=10: ",PorTip10YuNS4," %"),
     cex=1.3)

h <- hist(TipYulTS7,
          main="Tip Number from Yule simulations TreeSim\n( numSimb=500, age=6, lambda=1)",
          xlab="Tips")
text(x=quantile(h$breaks,0.7),
     y=quantile(h$counts,0.7),
     labels=paste("Tips>=10: ",PorTip10YuTS7," %"),
     cex=1.3)

h <- hist(TipYulNS7,
          main="Tip Number from Yule simulations NetSim\n( numSimb=500, age=6, lambda=1)",
          xlab="Tips")
text(x=quantile(h$breaks,0.7),
     y=quantile(h$counts,0.7),
     labels=paste("Tips>=10: ",PorTip10YuNS7," %"),
     cex=1.3)

#dev.off()
