arg1<-commandArgs(trailingOnly=TRUE)[1]
arg2<-commandArgs(trailingOnly=TRUE)[2]
arg3<-commandArgs(trailingOnly=TRUE)[3]
arg4<-as.numeric(commandArgs(trailingOnly=TRUE)[4])
eval(parse(text=arg1))
eval(parse(text=arg2))

#pat<-getwd()
#pathR
setwd(arg3)

library(SiPhyNetwork)#library to simulate Networks
#library(ape)#its a SiPhyNetwork dependence, and it is activated when you activate SiPhyNetwork or geiger
library(geiger)#is.phylo()##Information about tree goes extinct=0 and no extinct tips are sampled=1

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
if(i==2){res<-"Root"}else if(i==3){res<-"NoRoot"}else if (i<=1){res<-"No-Identificable"}
return(res)
}

#This function use the functions described before and determinate if the sibling reticulation can be Recognozible by snaq
Identiff<-function(Tree2,Hib1,Hib2,ED,father,RootTest,allret){
Tiprev<-Tip(ED,reff=RootTest,father=father,allret=allret)
Roo<-AssRoot(ED,reff=RootTest,allret=allret)
if(any(Roo%in%c("No-Identificable",NA,NULL))){resAll<-"No-Identificable"}else{

NodTip1<-Tiprev=="NoTip"|Roo=="NoRoot"
NodTip2<-Tip(ED,reff=Hib1,father=father,allret=allret)=="NoTip"
NodTip3<-Tip(ED,reff=Hib2,father=father,allret=allret)=="NoTip"
res1<-sum(c(NodTip1,NodTip2,NodTip3),na.rm = TRUE)#using (+) whit NA or NULL error will be apear, but if we use sum(,na.rm = TRUE), we omith that.7
#res1<-NodTip1+NodTip2+NodTip3
if(res1>=2){resAll<-"Identificable"}else {resAll<-"No-Identificable"} }
return(resAll)
}


# This function analize the tree, determinate if there are reticulations and in the case that it exists evalaute if the reticulation is among siblings, and in the case that it exists evaluate if it is or not recognozible by snaq.

#Mod1
SibCross<-function(Tree2){
#Tree2<-bdNS[[1]]
EDche<-Tree2$"edge"
RetRvCiclChe<-Tree2$"reticulation"
RetRvCiclChe1<-RetRvCiclChe
##


if(nrow(RetRvCiclChe)==0){res<-list(infor="NotRet",ret=NA)}else{


#Buscar si las retocilaciones se muestran en el file de los nodos, si es asi, es una reticulacion ciclica
out2<-c()
for(i in 1:nrow(RetRvCiclChe)){
Retche<-RetRvCiclChe[i,]
out<-c()
for(k in 1:nrow(EDche)){
re<-all(Retche==EDche[k,])
out<-c(out,re)
}
res<-any(out)
out2<-c(out2,res)
}

#EDcheWioutCicli
#Lista de reticulaciones ciclicas
RetClicli<-matrix(RetRvCiclChe[out2,],ncol=2)

#Remover las reticulaciones ciclicas del file de nodos y de las reticulaciones 
if(nrow(RetClicli)==0){
RetEv<-RetRvCiclChe
ED<-EDche
}else{
ED<-EDche
RetEv<-matrix(RetRvCiclChe[!out2,],ncol=2)
}


RetEv<-RetEv
ED<-ED
allret<-unique(as.numeric(RetEv))


##This function jumps until it reaches the nearest node, and in the way record the reticulation nodes which are no level one
out3<-c()
for(i in 1:length(allret)){
reff=allret[i]

ti12<-reff
Log=TRUE
out<-c()
while(Log){
ti12<-ED[ED[,1]==ti12,2][1]
Log<-ti12%in%allret
ro<-allret[allret%in%ti12]
out<-c(out,ro)
}
if(length(out)==0){res<-NA}else{res<-c(reff,out)}
out3<-c(out3,res)
}



NotL1<-na.omit(unique(out3))
if(length(NotL1)==0){RetEv<-RetEv;RetNoLev1=matrix(ncol=2)[-1,]}else{
logi1<-!(RetEv[,1]%in%NotL1|RetEv[,2]%in%NotL1)
RetNoLev1<-matrix(RetEv[!logi1,],ncol=2)
RetEv<-matrix(RetEv[logi1,],ncol=2)}


###############################################
ED<-ED
RetEv<-RetEv
allret<-unique(as.numeric(RetEv))
#######







if(nrow(RetEv)==0){
res<-list(infor=c(rep("not level-1",nrow(RetNoLev1)),rep("RetCiclicI",nrow(RetClicli))),ret=rbind(RetClicli,RetNoLev1))}else{
outRet<-c()

for(h in 1:nrow(RetEv)){
RetEv1<-RetEv[h,]
Hib1<-RetEv1[2]#Hibrid
Hib2<-RetEv1[1]#Father
ED<-Tree2$edge

Or1whitoutRet<-AssJupRet(ED,Hib1,allret)#Find the father omiting the reticulations
Or2whitoutRet<-AssJupRet(ED,Hib2,allret)#Find the father omiting the reticulations

Orig1<-ED[ED[,2]==Hib1,1]#Hibrid
Orig2<-ED[ED[,2]==Hib2,1]#Father
 
#res<-Orig1==Orig2

if(Or1whitoutRet==Or2whitoutRet){
father<-Orig1
RootTest<-AssJupRet(ED,Orig1,allret)#Go direction to the root avoiding the reticulation nodes.

Ident<-Identiff(Tree2,Hib1=Hib1,Hib2=Hib2,ED,father,RootTest,allret=allret)
SiblingCross<-Ident
}else{SiblingCross<-"Identificable"}#"NoSibCross"

outRet<-c(outRet,SiblingCross)
}
outRet<-c(outRet,c(rep("not level-1",nrow(RetNoLev1)),rep("RetCiclicI",nrow(RetClicli))))

RetEv<-rbind(RetEv,rbind(RetClicli,RetNoLev1))
res<-list(infor=outRet,ret=RetEv)
}


}




return(res)
}
######
RetAnalysis<-function(TreeAnalize1){
mat<-matrix(0,ncol=4,nrow=1)
colnames(mat)<-c("NotIdent", "Ident", "NotLevel-1","RetCiclicI")
colnam<-c("No-Identificable","Identificable","not level-1","RetCiclicI")

#plot(TreeAnalize1)
sibcro1<-SibCross(TreeAnalize1)
tab1<-table(sibcro1$infor)
namtab<-names(tab1)
for(j in 1:length(namtab)){
mat[,colnam==namtab[j]]<-tab1[j]
}
return(mat)
}






NetDescription<-function(bdNS){
out1<-c()
for(k in 1:length(bdNS)){
Tree2<-bdNS[[k]]
###Extracting general Information
meEdgeL<-mean(Tree2$"edge.length")#mean of edge length
sdEdgeL<-var(Tree2$"edge.length")
cvEdgeL<-sdEdgeL/meEdgeL#cv of edge length
#TipNum<-length(Tree2$"tip.label")# Number of tips
#PorExtin<-(TipNum-10)/TipNum#Proportions of extintct taxa 
NodeNum<-Tree2$Nnode#Number of internal Nodes# Maybe we can also divide it for TipNum
LT2<-ltt(Tree2,plot=FALSE)
Tim<-LT2$times
TimeTree<-max(Tim)#Time to reach 10 taxa
#extracting information about the Reticulations
Ret<-Tree2$"reticulation"
Ret1<-Ret
RetNum<-nrow(Ret)# Reticulation Number
if(RetNum>=1){
for(i in 1:nrow(Ret1)){
for(j in 1:ncol(Ret1)){
Ret1[i,j]<-Tim[names(Tim)%in%Ret1[i,j]]
}
}
TimeStartHibr<-apply(Ret1,1,max)/TimeTree#Time when birth the hibirds, we can divide it for the time to reach 10 taxa (TimeTree)
TSHmean<-mean(TimeStartHibr)# Average time of birth of hybrids 
TSHmin<-min(TimeStartHibr)# minimum time of birth of hybrids 
TSHmax<-max(TimeStartHibr)# maximun time of birth of hybrids 
}else{
TSHmean<-0
TSHmin<-0
TSHmax<-0
}

SiblingCross<-RetAnalysis(TreeAnalize1=Tree2)


#If reticulation number is one, TSHmean, TSHmin and TSHmax are the same, if there are not reticulation the value is cero.
#dat<-data.frame(TipNum,NodeNum,PorExtin,meEdgeL,cvEdgeL,TimeTree,RetNum,TSHmean,TSHmin,TSHmax)
dat<-data.frame(NodeNum,meEdgeL,cvEdgeL,TimeTree,RetNum,TSHmean,TSHmin,TSHmax,SiblingCross)
out1<-rbind(out1,dat)
}
return(out1)
}
###
#hybprops Combinations
Upp1<-1
Low1<-0
Gen<-c(rep(Upp1,4),rep(Low1,4))
Deg<-rep(c(Upp1,Upp1,Low1,Low1),2)
Neu<-rep(c(Low1,Upp1),4)


HybTab1<-data.frame(Gen,Deg,Neu)
HybTab<-HybTab1[-7,]#No acepta la combinacion c(0,0,0)


###

n1<-n1

numbsim1<-numbsim1
set.seed(arg4)
bdNS <- sim.bdh.taxa.ssa(n=n1, 
                                 numbsim = numbsim1, 
                                 lambda = lambda, 
                                 mu = mu , 
                                 nu = nu,
                                 hybprops = hybprops, 
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
bdNS <- bdNS[sapply(X = bdNS, FUN = is.phylo)]

out<-c()
for(i in 1:length(bdNS)){
res<-SibCross(Tree2=bdNS[[i]])$infor
ro<-any(res%in%"not level-1")
out<-c(out,ro)
}


id<-which(out)

for(i in 1:length(id)){
tree<-bdNS[[id[i]]]
res<-SibCross(Tree2=tree)
#plot(tree)
#nod<-res$ret[res$infor=="not level-1",]
#nodelabels(node=nod)
write.net(tree,file=paste("RNetwork_",i,sep=""))
}
