### R functions needed for generate_datasets.R
### GAB, CA, CSL (August 2022)


## GAB: original construction of the is.phylo function in geiger. 
## Unnecessary to load the whole thing for just a function
is.phylo <- function(x) {
  "phylo" %in% class(x)
}

## GAB:: a version of is.rooted that works on networks
# phy: an object that inherits from the class phylo.
# returns TRUE if the only node which is an ancestor is not a descendant itself and it
# has exactly two children 
is.really.rooted <- function(phy) {
    root_candidates <- phy$edge[!(phy$edge[,1] %in% phy$edge[,2]), 1]
    anc_table <- table(root_candidates)
    if (min(anc_table) == 2) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

library(SiPhyNetwork) #library to simulate Networks, this use "ape" as dependence
library(ape)

set.seed(6550)
networks <- sim.bdh.taxa.ssa(n = 15,
                             numbsim = 150,
                             lambda = 0.9,
                             mu = 0 ,
                             nu = 0.4,#0.02
                             hybprops = c(1, 1, 1),
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

net1 <- networks[[1]]

plot(net1)
nodelabels(cex=0.7) #add node numbers
tiplabels(cex=0.7) 






# commenting scripts ####

# sibling_ret ####

# input: Evolutionary network
# Output: message "No sibling reticulation" or "Sibling reticulation"
# and also shows the No sibling reticulation nodes.

# Find if there are reticulation among siblings, that is not recognozible by snaq
# it happens when the father of the nodes of the reticulations are a root
# or when the son of both reticulation node are a tip

sibling_ret <- function(network){
  ## Estimations ####
  ### rets and edges ####
  rets <- network$reticulation# Reticulations
  ed_1 <- network$edge# Edges
  
  ### Tips ####
  # A tip is a node that has no children
  tips <- ed_1[,2][!ed_1[,2]%in%ed_1[,1]]
  
  ### root ####
  # padre que no es hijo
  p1 <- ed_1[!ed_1[,1]%in%ed_1[,2],1]
  h_p1 <- ed_1[ed_1[,1]%in%p1,2]# hijos  de p1
  # si solo tiene un hijo ese hijo seria el verdadero root
  # si tiene dos hijo el sería el root
  if(length(h_p1)==1){root <- h_p1}
  if(length(h_p1)==2){root <- p1}
  root
  
  ## body ####
  out <- c()
  for(i in 1:nrow(rets)){
    #i=1
    ret_i <- rets[i,]
    
    # parents of reticulations nodes
    # parent of the start node of the reticulation
    father_ret1 <- ed_1[ed_1[,2]==ret_i[1],1]
    # parent of the start node of the reticulation
    father_ret2 <- ed_1[ed_1[,2]==ret_i[2],1] 
    
    # if the reticulation node have the same faher->sibling cross
    sibling_log1 <-father_ret1==father_ret2
    
    # The childs of the reticulation node are tips?
    child.ret.node_are.tips <- all(ed_1[,2][ed_1[,1]%in%ret_i]%in%tips)
    
    # patens of both reticulation nodes are root?
    parents.ret.node_are.root <- father_ret1==root
    
    if((sibling_log1&child.ret.node_are.tips)|(sibling_log1&parents.ret.node_are.root)){
      out <- rbind(out,ret_i)
    }else{}
    
  }
  
  if(length(out)==0){
    res<-list(mess= "No sibling reticulation", rets="")
  }else{
    res<-list(mess= "Sibling reticulation", rets=out)
  }
  return(res)
}





# ret_higher_level_detection complementary functions ####
# This function finds the start node of the edge that contains the reticulation
# input: Reticulation node to find the limits of the edge and the network
# Output: Start node of the edge that contains the reticulation
start_node_edge_ret <- function(ret_to_evaluate,network){
  edges<-network$edge
  rets<-network$reticulation
  
  logi_father_is_ret <- TRUE
  rets_all <- as.numeric(rets)
  while(logi_father_is_ret){
    father_ret1 <- edges[edges[,2]==ret_to_evaluate, 1]
    logi_father_is_ret <- any(father_ret1==rets_all)
    logi_father_is_ret
    ret_to_evaluate<-father_ret1
  }
  father_node_rets<-father_ret1
  
  return(father_node_rets)
}

# This function finds the end node of the edge that contains the reticulation
# input: Reticulation node to find the limits of the edge and the network
# Output: End node or leave of the edge that contains the reticulation
end_node_edge_ret <- function(ret_to_evaluate,network){
  edges<-network$edge
  rets<-network$reticulation
  leaves <- edges[,2][!edges[,2]%in%edges[,1]]
  
  logi_son_is_ret <- TRUE
  rets_all <- as.numeric(rets)
  if(any(ret_to_evaluate==leaves)){son_node_rets <- ret_to_evaluate}else{
    
    while(logi_son_is_ret){
      son_ret1 <- edges[edges[,1]==ret_to_evaluate, 2]
      logi_son_is_ret <- !(!any(son_ret1==rets_all)|any(son_ret1%in%leaves))
      logi_son_is_ret
      ret_to_evaluate <- son_ret1
    }
    son_node_rets <- son_ret1
  }
  son_node_rets
  
  return(son_node_rets)
}



# Detection of higher level reticulation ####
# input: Evolutionary network
# output: message of "level_one" or "higher_level"
# if there are at least one higher_level reticulation
# show also the reticulations and the edges with the number of reticulations 
# so we can see th level of reticulations

ret_higher_level_detection <- function(network){
  ## Estimations ####
  ### Reticualtions and edges
  net1 <- network
  rets <- net1$reticulation
  edges <- net1$edge
  
  ### sons and fathers of the reticulations ####
  hijos <- edges[edges[,1]%in%rets,2]
  padres <- edges[edges[,2]%in%rets,1]
  hijos_padres_ret <- c(hijos, padres)
  
  # If some of the child or parent nodes of the reticulations belong to the list 
  # of total nodes in the edges, it means that there is at least one branch with 
  # more than one reticulation (higher level reticulation)
  not_single_log <- any(hijos_padres_ret %in% rets)# Logic
  
  branch_rets <-c()
  if(not_single_log){
    
    
    nodes_ret_higher_level <- hijos_padres_ret[hijos_padres_ret %in% rets]
    rets_hig_log <- rets[,1]%in%nodes_ret_higher_level|rets[,2]%in%nodes_ret_higher_level
    rets_hig <- rets[rets_hig_log,]

    for(i in 1:length(nodes_ret_higher_level)){
      # Find the limits of the branch that contains the reticulation nodes 
      # so we can count the reticulation nodes inside the branch and 
      # determiante the level 
      ret_to_find_limits <- nodes_ret_higher_level[i]
      branch_ini <- start_node_edge_ret(ret_to_find_limits , network=net1)
      branch_fini <- end_node_edge_ret(ret_to_find_limits, network=net1)
      res <- data.frame(Node=ret_to_find_limits, branch_ini = branch_ini, branch_fini = branch_fini)
      branch_rets <- rbind(branch_rets,res);branch_rets 
    }
    
    # 
    branch_rets$Num_rets<-1
    banch_freq_rets <- aggregate(Num_rets~branch_ini+branch_fini,data=branch_rets,FUN="length")
    
    list1<-list(resp ="higher_level", rets=rets_hig, branch =banch_freq_rets)
  }else{
    list1<-list(resp ="level_one", rets="",branch="")
  }
  
  return(list1)
}



#is_not_single_ret(networks[[1]])

#sibling_ret(networks[[1]])






###########################################################################################
#################################################################################################
##################################################################################################

#We want to know if the reticulation node produce a Tip or Node
##This function goes in the direction of the tips omiting reticulations nodes until it finds one no reticulation node
desJupRet <- function(ED, reff, father, allret){
  ti12 <- reff
  Log <- TRUE
  while(Log){
    ti12 <- ED[ED[, 1] == ti12 & ED[, 2]!=father, 2][1]
    #we use father because when we want to evaluate RootRet we want to know if [Tip or Node 3] is a node or tip, no go for the reticulation way.
    Log <- ti12 %in% allret
  }
  return(ti12)
}


#This function goes in the direction of the tips jumping node by node and determinate if the refference (reff) produce a Tip or not
tip <- function(ED, reff, father, allret){
  Log <- TRUE
  i <- 0
  while(Log){
    reff <- desJupRet(ED=ED, reff=reff, father=father, allret=allret)
    #
    Log <- !(i == 2|is.na(reff))
    i <- i+1
  }
  if(i == 2){
      res <- "Tip"
  } else if(i == 3){
      res <- "NoTip"
  }
  return(res)
}


##This function goes in the direction of the root omiting reticulations nodes until it finds one no reticulation node
assJupRet <- function(ED, reff, allret){
  ti12 <- reff
  Log=TRUE
  while(Log){
    ti12 <- ED[ED[, 2] == ti12, 1][1]###Loop
    Log <- ti12 %in% allret
  }
  return(ti12)
}

#This function goes in the direction of the root jumping node by node and determinate if the refference (reff) produce a root or node
assRoot <- function(ED, reff, allret){
  Log <- TRUE
  i <- 0
  while(Log){
    reff <- assJupRet(ED=ED, reff=reff, allret=allret)
    Log <- !(i == 2|is.na(reff))
    i <- i+1
  }
  if(i == 2){
      res <- "Root"
  } else if(i == 3){
      res <- "NoRoot"
  } else if (i<=1){
      res <- "No-Identificable"
  }
  return(res)
}

#This function use the functions described before and determinate if the sibling reticulation can be Recognozible by snaq
identiff <- function(Tree2, Hib1, Hib2, ED, father, RootTest, allret){
  Tiprev <- tip(ED, reff=RootTest, father=father, allret=allret)
  Roo <- assRoot(ED, reff=RootTest, allret=allret)
  if(any(Roo %in% c("No-Identificable", NA, NULL))){
      resAll <- "No-Identificable"
  } else {
    NodTip1 <- Tiprev == "NoTip" | Roo == "NoRoot"
    NodTip2 <- Tip(ED, reff=Hib1, father=father, allret=allret) == "NoTip"
    NodTip3 <- tip(ED, reff=Hib2, father=father, allret=allret) == "NoTip"
    res1 <- sum(c(NodTip1, NodTip2, NodTip3), na.rm = TRUE)#using (+) whit NA or NULL error will be apear, but if we use sum(, na.rm = TRUE), we omith that.7
    #res1 <- NodTip1+NodTip2+NodTip3
    if(res1 >= 2){
        resAll <- "Identificable"
    } else {
        resAll <- "No-Identificable"
    }
  }
  return(resAll)
}


# This function analize the tree, determinate if there are reticulations and in the case that it exists evalaute if the reticulation is among siblings, and in the case that it exists evaluate if it is or not recognozible by snaq.

#Mod1
sibCross <- function(Tree2){
  #Tree2 <- bdNS[[1]]
  EDche <- Tree2$"edge"
  RetRvCiclChe <- Tree2$"reticulation"
  RetRvCiclChe1 <- RetRvCiclChe

    if(nrow(RetRvCiclChe) == 0){
        res <- list(infor = "NotRet", ret = NA)
    } else { #GAB esto no cierra?
    #Buscar si las reticulaciones se muestran en el file de los nodos, si es asi, es una reticulacion ciclica
    out2 <- c()
    for(i in 1:nrow(RetRvCiclChe)){
      Retche <- RetRvCiclChe[i, ]
      out <- c()
      for(k in 1:nrow(EDche)){
        re <- all(Retche == EDche[k, ])
        out <- c(out, re)
      }
      res <- any(out)
      out2 <- c(out2, res)
    }

    #EDcheWioutCicli
    #Lista de reticulaciones ciclicas
    RetClicli <- matrix(RetRvCiclChe[out2, ], ncol=2)

    #Remover las reticulaciones ciclicas del file de nodos y de las reticulaciones
    if(nrow(RetClicli) == 0){
      RetEv <- RetRvCiclChe
      ED <- EDche
    } else {
      ED <- EDche
      RetEv <- matrix(RetRvCiclChe[!out2, ], ncol=2)
    }


    RetEv <- RetEv
    ED <- ED
    allret <- unique(as.numeric(RetEv))


    ##This function jumps until it reaches the nearest node, and in the way record the reticulation nodes which are no level one
    out3 <- c()
    for(i in 1:length(allret)){
      reff=allret[i]
      ti12 <- reff
      Log <- TRUE
      out <- c()
      while(Log){
        ti12 <- ED[ED[, 1] == ti12, 2][1]
        Log <- ti12%in%allret
        ro <- allret[allret%in%ti12]
        out <- c(out, ro)
      }
      if(length(out) == 0){res <- NA}else{res <- c(reff, out)}
      out3 <- c(out3, res)
    }

    NotL1 <- na.omit(unique(out3))
    if(length(NotL1) == 0){RetEv <- RetEv;RetNoLev1=matrix(ncol=2)[-1, ]}else{
      logi1 <- !(RetEv[, 1]%in%NotL1|RetEv[, 2]%in%NotL1)
      RetNoLev1 <- matrix(RetEv[!logi1, ], ncol=2)
      RetEv <- matrix(RetEv[logi1, ], ncol=2)}


    ###############################################
    #GAB qué hacen las siguientes dos líneas de código?
    ED <- ED
    RetEv <- RetEv
    allret <- unique(as.numeric(RetEv))
    #######


    if(nrow(RetEv) == 0){
      res <- list(infor=c(rep("not level-1", nrow(RetNoLev1)), rep("RetCiclicI", nrow(RetClicli))), ret=rbind(RetClicli, RetNoLev1))}else {
        outRet <- c()

        for(h in 1:nrow(RetEv)){
          RetEv1 <- RetEv[h, ]
          Hib1 <- RetEv1[2]#Hibrid
          Hib2 <- RetEv1[1]#Father
          ED <- Tree2$edge

          Or1whitoutRet <- assJupRet(ED, Hib1, allret)#Find the father omiting the reticulations
          Or2whitoutRet <- assJupRet(ED, Hib2, allret)#Find the father omiting the reticulations

          Orig1 <- ED[ED[, 2] == Hib1, 1]#Hibrid
          Orig2 <- ED[ED[, 2] == Hib2, 1]#Father

          #res <- Orig1 == Orig2

          if(Or1whitoutRet == Or2whitoutRet){
            father <- Orig1
            RootTest <- assJupRet(ED, Orig1, allret)#Go direction to the root avoiding the reticulation nodes.

            Ident <- identiff(Tree2, Hib1=Hib1, Hib2=Hib2, ED, father, RootTest, allret=allret)
            SiblingCross <- Ident
          }else{SiblingCross <- "Identificable"}#"NoSibCross"

          outRet <- c(outRet, SiblingCross)
        }
        outRet <- c(outRet, c(rep("not level-1", nrow(RetNoLev1)), rep("RetCiclicI", nrow(RetClicli))))

        RetEv <- rbind(RetEv, rbind(RetClicli, RetNoLev1))
        res <- list(infor=outRet, ret=RetEv)
      }


  }

  return(res)
}
