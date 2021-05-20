library(NetSim) # simulate bd time-tree networks: https://github.com/jjustison/NetSim
library(TreeSim) # simulate bd time-trees
library(ape) # tools for handling phylo objects
library(diversitree) # tools for estimating lambda and mu
library(phytools)
library("dplyr")    # Data manipulation
library(ggplot2)
library(corrplot)
library(cluster)

library(NetSim) # simulate bd time-tree networks: https://github.com/jjustison/NetSim
library(geiger)
library(ape) # tools for handling phylo objects



set.seed(200)
numbsim=400
n1<-6
bdNS <- NetSim::sim.bdh.taxa.ssa(n=n1, 
                                 numbsim = numbsim, 
                                 lambda = 0.9, 
                                 mu = 0, 
                                 nu = 0.3,
                                 hybprops = c(0.5,0.5,0.5), 
                                 hyb.inher.fxn = NetSim::make.beta.draw(1, 1), 
                                 frac = 1,
                                 mrca = FALSE, 
                                 complete = TRUE, 
                                 stochsampling = FALSE, 
                                 hyb.rate.fxn = NULL, 
                                 trait.model = NULL)
#
bdNS <- bdNS[!sapply(X = bdNS, FUN = is.null)]
#Information about tree goes extinct=0 and no extinct tips are sampled=1
bdNS <- bdNS[sapply(X = bdNS, FUN = is.phylo)]#birth dead tree with NetSim with value 0 or 1
 


#Squeme of evaluation if a sibling Reticulation is Recognosible
 

#			
#				         _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  _ _ _ _ _ _ _ _[Tip or Node 4]     
#				        |
#				        |
#				        |
#                                             Root_ _ _ _[Node#]               _ _ _ _ _ _ _ _ _ _ _ _ _ _ _[Tip or Node 3]
#				        |		|
#                                                           |                 |
#				        |		|
#				        |_ _ _ _ _ _ _RootRet			   
#						|        _ _ _ RetNode1_ _ _ _[Tip or Node 2]
#						|       |         |
#						|_ _ _father      |
#						         |_ _ _ RetNode2_ _ _ _[Tip or Node 1]	

#I evaluate if a RootRet produce a Tip or Node			
## I considered that RootRet produce a "Node" if [Tip or Node 3] is a node, or [Node#] is a node.
## On the other hands if [Node#] is a root and [Tip or Node 3] is a tip, I considered that RootRet produce a "Tip"

#Then I evaluate if a RetNode1 and RetNode2 produce a Tip or Node

#Finaly if number of tips are >=2, the siblign reticulation is not Recognosible otherwise is Recognosible.



#We want to know if the reticulation node produce a Tip or Node
##This function goes in the direction of the tips omiting reticulations nodes until it finds one no reticulation node
DesJupRet<-function(ED,reff,father,allret){
ti12<-reff
Log=TRUE
while(Log){
ti12<-ED[ED[,1]==ti12&ED[,2]!=father,2][1]
#we use father because when we want to evaluate RootRet we want to know if [Tip or Node 3] is a node or tip, no go for the reticulation way.
Log<-ti12%in%allret
}
return(ti12)
}


#This function goes in the direction of the tips jumping node by node and determinate if the refference (reff) produce a Tip or not
Tip<-function(ED,reff,father,allret){
Log=TRUE
i=0
while(Log){
reff<-DesJupRet(ED=ED,reff=reff,father=father,allret=allret)
#
Log<-!(i==2|is.na(reff))
i=i+1
}
if(i==2){res<-"Tip"}else if(i==3){res<-"NoTip"}
return(res)
}


##This function goes in the direction of the root omiting reticulations nodes until it finds one no reticulation node
AssJupRet<-function(ED,reff,allret){
ti12<-reff
Log=TRUE
while(Log){
ti12<-ED[ED[,2]==ti12,1][1]###Loop
Log<-ti12%in%allret
}
return(ti12)
}

#This function goes in the direction of the root jumping node by node and determinate if the refference (reff) produce a root or node
AssRoot<-function(ED,reff,allret){
Log=TRUE
i=0
while(Log){
reff<-AssJupRet(ED=ED,reff=reff,allret=allret)
Log<-!(i==2|is.na(reff))
i=i+1
}
if(i==2){res<-"Root"}else if(i==3){res<-"NoRoot"}#else if (i<=1){res<-"No-Identificable"}
return(res)
}

#This function use the functions described before and determinate if the sibling reticulation can be Recognozible by snaq
Identiff<-function(Tree2,Hib1,Hib2,ED,father,RootTest,allret){
Tiprev<-Tip(ED,reff=RootTest,father=father,allret=allret)
Roo<-AssRoot(ED,reff=RootTest,allret=allret)
#if(any(Roo%in%c("No-Identificable",NA,NULL))){resAll<-"No-Identificable"}else{

NodTip1<-Tiprev=="NoTip"|Roo=="NoRoot"
NodTip2<-Tip(ED,reff=Hib1,father=father,allret=allret)=="NoTip"
NodTip3<-Tip(ED,reff=Hib2,father=father,allret=allret)=="NoTip"
#res1<-sum(c(NodTip1,NodTip2,NodTip3),na.rm = TRUE)#using (+) whit NA or NULL error will be apear, but if we use sum(,na.rm = TRUE), we omith that.7
res1<-NodTip1+NodTip2+NodTip3
if(res1>=2){resAll<-"Identificable"}else {resAll<-"No-Identificable"} #}
return(resAll)
}


# This function analize the tree, determinate if there are reticulations and in the case that it exists evalaute if the reticulation is among siblings, and in the case that it exists evaluate if it is or not recognozible by snaq.

SibCross<-function(Tree2){
RetEv<-Tree2$"reticulation"
allret<-unique(as.numeric(RetEv))
ED<-Tree2$edge

if(nrow(RetEv)==0){SiblingCross<-"NotRet"
res<-list(infor=SiblingCross,ret=NA)}else{
outRet<-c()

for(h in 1:nrow(RetEv)){
RetEv1<-RetEv[h,]
Hib1<-RetEv1[2]#Hibrid
Hib2<-RetEv1[1]#Father
ED<-Tree2$edge
Orig1<-ED[ED[,2]==Hib1,1]#Hibrid
Orig2<-ED[ED[,2]==Hib2,1]	#Father 
res<-Orig1==Orig2
if(res){
father<-Orig1
RootTest<-ED[ED[,2]==Orig1,1]

Ident<-Identiff(Tree2,Hib1=Hib1,Hib2=Hib2,ED,father,RootTest,allret=allret)
SiblingCross<-Ident
}else{SiblingCross<-"NoSibCross"}

outRet<-c(outRet,SiblingCross)
}
res<-list(infor=outRet,ret=RetEv)
}

return(res)
}


#Tree2<-bdNS[[13]]
#SibCross(Tree2)


#No Identificable for snaq
out<-c()
for(i in 1:length(bdNS)){
res<-SibCross(bdNS[[i]])$infor
res2<-any(res=="Identificable")
out<-c(out,res2)
}


posNotIdent1<-which(out);posNotIdent1
TreeAnalize1<-bdNS[[posNotIdent1[4]]]

sibcro1<-SibCross(TreeAnalize1)
nod1<-sibcro1$ret[sibcro1$infor=="Identificable",]
plot(TreeAnalize1)
nodelabels(node=nod1)

#plot(plottable.net(TreeAnalize1))
#nodelabels(node=nod1)





#Identificable for snaq
out2<-c()
for(i in 1:length(bdNS)){
res<-SibCross(bdNS[[i]])$infor
res2<-any(res=="No-Identificable")
out2<-c(out2,res2)
}
posNotIdent<-which(out2);posNotIdent

TreeAnalize<-bdNS[[posNotIdent[2]]]

sibcro<-SibCross(TreeAnalize)
nod<-sibcro$ret[sibcro$infor=="No-Identificable",]
plot(TreeAnalize)
nodelabels(node=nod)


#plot(plottable.net(TreeAnalize))
#nodelabels(node=nod)





