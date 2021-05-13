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







DesJupRet<-function(ED,reff,father,allret){
ti12<-reff
Log=TRUE
while(Log){
ti12<-ED[ED[,1]==ti12&ED[,2]!=father,2][1]###Loop
Log<-ti12%in%allret
}
return(ti12)
}


Tip<-function(ED,reff,father,allret){
Log=TRUE
i=0
while(Log){
reff<-DesJupRet(ED=ED,reff=reff,father=father,allret=allret)
Log<-!(i==2|is.na(reff))
i=i+1
}
if(i==2){res<-"Tip"}else if(i==3){res<-"NoTip"}
return(res)
}


AssJupRet<-function(ED,reff,allret){
ti12<-reff
Log=TRUE
while(Log){
ti12<-ED[ED[,2]==ti12,1][1]###Loop
Log<-ti12%in%allret
}
return(ti12)
}

AssRoot<-function(ED,reff,allret){
Log=TRUE
i=0
while(Log){
reff<-AssJupRet(ED=ED,reff=reff,allret=allret)
Log<-!(i==2|is.na(reff))
i=i+1
}
if(i==2){res<-"Root"}else if(i==3){res<-"NoRoot"}
return(res)
}



Identiff<-function(Tree2,Hib1,Hib2,ED,father,RootTest,allret){
Tiprev<-Tip(ED,reff=RootTest,father=father,allret=allret)
Roo<-AssRoot(ED,reff=RootTest,allret=allret)
NodTip1<-Tiprev=="NoTip"|Roo=="NoRoot"
NodTip2<-Tip(ED,reff=Hib1,father=father,allret=allret)=="NoTip"
NodTip3<-Tip(ED,reff=Hib2,father=father,allret=allret)=="NoTip"
res1<-NodTip1+NodTip2+NodTip3
if(res1>=2){resAll<-"Identificable"}else {resAll<-"No-Identificable"}
resAll
}



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


#No Identificable for snaq
out<-c()
for(i in 1:length(bdNS)){
res<-SibCross(bdNS[[i]])$infor
res2<-any(res=="No-Identificable")
out<-c(out,res2)
}
which(out)

posNotIdent1<-which(out);posNotIdent1
TreeAnalize1<-bdNS[[posNotIdent1[2]]]

sibcro1<-SibCross(TreeAnalize1)
nod1<-sibcro1$ret[sibcro1$infor=="No-Identificable",]
plot(TreeAnalize1)
nodelabels(node=nod1)

#plot(plottable.net(TreeAnalize1))
#nodelabels(node=nod1)



#Identificable for snaq
out2<-c()
for(i in 1:length(bdNS)){
res<-SibCross(bdNS[[i]])$infor
res2<-any(res=="Identificable")
out2<-c(out2,res2)
}
posNotIdent<-which(out2);posNotIdent

TreeAnalize<-bdNS[[posNotIdent[4]]]

sibcro<-SibCross(TreeAnalize)
nod<-sibcro$ret[sibcro$infor=="Identificable",]
plot(TreeAnalize)
nodelabels(node=nod)


#plot(plottable.net(TreeAnalize))
#nodelabels(node=nod)





