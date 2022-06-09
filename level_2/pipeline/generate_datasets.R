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
















######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

####GAB Parameters were specified on notion in https://www.notion.so/Project-on-SNaQ-limitations-80e419d6d58a42308e4bf0f5cceb9752#e5df0d1a4c97433ca14f2e6388be966a
numbsim <- 150
#numbsim <- 50
ntips <- c(15, 30, 50)
#ntips <- c(15)
lambda <- 0.9
mu <- 0
nu <- c(0.02, 0.04)
#nu <- c(0.02)
hybprops <- c(1, 1, 1) # we are starting with all-equal
ngt <- c(100, 1000, 10000)
#ngt <- c(50)
set.seed(2022) # our initial seed, all the others below will depend on sample.int()

#GAB  original construction of the is.phylo function in geiger. Unnecessary to load the whole thing for just a function
is.phylo <- function(x) {
  "phylo" %in% class(x)
}

library(SiPhyNetwork) #library to simulate Networks, this use "ape" as dependence

#setwd(patscripts)
#GAB create a directory for storing the networks and their subdirs. the .. is necessary as it assumes that the script generate_datasets.R us run from the directory pipeline, which should only contain code and not data
dir.create(path = "../data")


#3.1 Simulate networks
#GAB set each seed as below, They are all pseudo-random, but dependent on the first set.seed(2022) above

#GAB nested for for generating the params one at a time for those which vary, i.e., ntips, nu, ngt

for (i in ntips) {
  for (j in nu) {
    for (k in ngt) {
      #i=ntips[1];j=nu[1];k=ngt[1]
      setwd("../data")
      r_seed <- sample.int(n = 1e6, size = 1)
      cat("r_seed = ", r_seed, "\n", sep = "")
      set.seed(r_seed)
      
      networks <- sim.bdh.taxa.ssa(n = i,
                                   numbsim = numbsim, 
                                   lambda = lambda, 
                                   mu = mu , 
                                   nu = j,
                                   hybprops = hybprops, 
                                   hyb.inher.fxn = make.beta.draw(1, 1), 
                                   frac = 1,
                                   mrca = FALSE, 
                                   complete = TRUE, 
                                   stochsampling = FALSE, 
                                   hyb.rate.fxn = NULL, 
                                   trait.model = NULL)
      #GAB Code for removing bad networks before writing
      # get rid of null trees which go extinct=0 and no extinct tips are sampled=1
      networks <- networks[!sapply(X = networks, FUN = is.null)]
      networks <- networks[sapply(X = networks, FUN = is.phylo)]
      
      #Code to select only networks (omit trees)
      file_networks <- vector(length = length(networks))
      for(x in 1:length(networks)){          
        file_networks[x] <- as.logical(nrow(networks[[x]]$reticulation))
      }
      networks <- networks[file_networks]
      
      net_counter <- 1
      for (y in networks) {
        
        #net_counter <- 3;y=networks[[net_counter]]
        filename <- paste("network", net_counter,
                          "_ntips_", i,
                          "_nu_", j,
                          "_ngt_", k, sep="")
        
        dir.create(filename)
        setwd(filename)
        
        # CA: Not level one warning ########
        res<-SibCross(Tree2=y)$infor
        ro<-any(res%in%"not level-1")
        if(ro){
          sink("not_level1_warn.txt")
          cat(res)
          sink()
        }
        ################################
        
        gt_seed <- sample.int(n = 1e6, size = 1)
        cat("gt_seed = ", gt_seed, "\n", sep = "")
        extnewick_filename <- paste(filename, ".extnewick", sep = "")
        
        # write to string file
        SiPhyNetwork::write.net(net = y, file = extnewick_filename)
        
        hybridlambda_filename <- paste(filename, ".hybridlambda", sep = "")
        
        # convert from extnewick to hybridlambda
        
        system(paste("julia ../../pipeline/extnewick2hybridlambda.jl ", 
                     extnewick_filename, " ", hybridlambda_filename, sep = ""))
        
        
        # run hybrid-Lambda on filename_extnewick and capture the output
        system(paste("hybrid-Lambda -spcu ", "'",
                     hybridlambda_filename,"'", " -num ", 
                     k," -seed ", gt_seed, " -o ",
                     hybridlambda_filename," > hybridlambda_output 2>&1",sep=""))
        
        notultram_bool <- sum(grepl(x = readLines("hybridlambda_output"),
                                    pattern = "ERROR: Non-ultrametric tree"))#for Gustavo's Hybrid lambda
        notultram_bool
        
        if(notultram_bool>0){#Is Not ultrametric
          #GAB the network file in extnewick will be called appending the parameter values as well as a counter for numbering each network from 1 no length(networks)
          setwd("../")
          system(paste("rm -rf ", filename, sep = ""))
        } else {
          file.remove("hybridlambda_output")
          setwd("../")
        }
        net_counter <- net_counter + 1
      }
    }
  }
}
