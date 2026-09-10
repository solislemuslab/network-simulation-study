library(phangorn)
library(stringr)
library(MCMCprecision)
library(MSCquartets)

#This is a function that takes each network from PhyloSketch, removes characters before the hybrid and replaces them with 
#LGT# followed by the parameters. 
top2metric = function(X){
  {
    netsplit = strsplit(X,'')[[1]]
    sep = which('#'==netsplit)
    for(i in 1:length(sep))
    {
      netsplit[ sep[i] - 1 ] = ''
      netsplit[ sep[i] + 1 ] = 'Z'
    }
    RemoveNetwork = paste0(netsplit,collapse = "")
  }
  Y = paste0(netsplit, collapse = "")
#Reading each network and counting edges/ hybrid nodes
  
  N = read.evonet(text = Y)
  loopsize = dim(N$edge)[1] + 1
  jk = length(unique(str_extract_all(Y, "#Z\\d+")[[1]]))
  
  lens = rgamma(loopsize + 2*jk, 1)
  hybrid = runif(jk, 0.2,0.8)
  hybrid = c(rbind(hybrid, 1 - hybrid))
  j=1;
  
  
  ###adding branch lengths
  for(i in LETTERS[1:loopsize])
  {
    Y = str_replace(Y, i, paste0(i, ':', lens[j]))
    j = j + 1
  } #This is adding a : for each edge
  for(k in 1:jk){
    Y = str_replace(Y,
                    paste0('\\)#Z', k),
                    paste0(')#LGT', k, ':', lens[j], '::', hybrid[2*k - 1]))
    j = j + 1
    
    Y = str_replace(Y, 
                    paste0('#Z', k),
                    paste0('#LGT', k, ':', lens[j], '::', hybrid[2*k]))
    j = j + 1
  } #adding :: and replacing Z for LGT# and the number of hybrid associated to it
    Y
}
##Result should print (((((((C:3)#LGT3:7::0.2,D:1)L:3,E:2)M:1)#LGT2:5::0.3)#LGT1:4::0.1,((#LGT3:9::0.8,B:2)I:3,A:1)H:2)G:1,(#LGT1:3::0.9,(#LGT2:6::0.7,F:3)P:1)O:3)r;
## as an example. Below is such a case. 

#Example
X ="(((((((E)N#H3,F)Q)T#H2)O#H1,(((N#H3,D)L,C)K,B)J)I,A)M,(H,(O#H1,(T#H2,G)S)R)P)r;"
Y = top2metric(X)
Y

#Just to clean the code up
path ='C:/Users/rober/OneDrive/Desktop/Level - 3 networks Newick Notation/'
path2master = 'C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/Level-3 Master File/masterfilelevel35000.jl'
path2masterSNaQ = 'C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/Level-3 Master File/SNaQMasterlevel35000.jl'


####Note this loop is for all networks in the semi-directed case
#######Sampling with the normal distribution from a folder network

folder = list.files(paste0(path,"3a"))
for(j  in 1:length(folder)){
  file = folder[j]
  X = readLines(paste0(path, "3a/",file))
  X = paste0(X, collapse ="" )
  YY = top2metric(X)
  
  ####Sampling gene trees from the given folder
  n = 1000 ### number of sample gene trees
  MS = readLines(paste0(path2master)) #Will read the master file in R and run this in Julia
  MS[9] = paste0('num_gts =' ,n, ';') #Create gene trees and names the file for n sample of gene trees
  MS[13] = paste0('newick =', "\"",YY,"\"") #Will paste the Newick for each gene tree
  MS[11] = paste0('graph_name =', "\"",file,"\"") #Will create the name after having the sample of n size
  MS[7] = paste0("cd(","\"","C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/Thesis_Julia_Project_Level3_5000/","3agenetrees","\"",")")
  write.table( MS ,paste0("C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/Juliafolder_Level3_5000/Generator3a/",file,".jl"),quote=F, append = FALSE, sep = " ", dec = ".", row.names = F, col.names = F) #Will format this in text for it to be transferred back to SNaQ. 
  system(paste0('julia ',paste0("C:\\Users\\rober\\OneDrive\\Desktop\\R_Project_Folder_Phylogenetics\\Juliafolder_Level3_5000\\Generator3a\\",file,".jl"))) #Run the command in Julia
  
  
}


####SNaQ
for(j  in 1:length(folder)){
  file = paste0("C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/Thesis_Julia_Project_Level3_5000/3agenetrees/",'gtrees_', folder[j])
  MS = readLines(paste0(path2masterSNaQ)) #Will read the master file in R and run this in Julia
  MS[2] = paste0('gt_location = string(',"\"",file,"\"",')')
  MS[3] = paste0('cd(',"\"",'C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/SNaQnetwork_Level3_5000/Generator3a',"\"",')')
  MS[15] = paste0('net1 = snaq!(starttree,CF, hmax=1, filename=',"\"",folder[j],'snaqrun',"\"",', seed=1234)')
  write.table( MS ,paste0("C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/Juliafolder_Level3_5000/Generator3a/",folder[j],'snaq',".jl"),quote=F, append = FALSE, sep = " ", dec = ".", row.names = F, col.names = F) #Will the format this in text for it to be transferred back to SNaQ. 
  system(paste0('julia ',paste0("C:\\Users\\rober\\OneDrive\\Desktop\\R_Project_Folder_Phylogenetics\\Juliafolder_Level3_5000\\Generator3a\\",folder[j],'snaq',".jl"))) #Run this command in Julia
  
}    

####Read Files SNaQ

for(j  in 1:length(folder)){
  fileSNAQ =  paste0("C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/SNaQnetwork_Level3_5000/Generator3a/",folder[j],'snaqrun.networks')
  
  #Goal of this is to read the text file, removes all the wording except to leave it in Newick Notation
  lineSNAQ = readLines(fileSNAQ)[1] #Reading the first line in the networks file
  lineSNAQ = gsub("\\, .*","",lineSNAQ) #Removes all wording and leaves it in Newick Notation 
  lineSNAQ = gsub(":::",":1:",lineSNAQ) #Removes the triple colon and replaces it with only 1 colon followed by the branch length of 1
  networkSNAQ = read.evonet(text=lineSNAQ)
  networkSNAQ$edge.length =  rep(1, length(networkSNAQ$edge.length)) #Adding edge lengths to each remaining taxa, length 1
  plot(networkSNAQ, main = paste0("network SNaQ", j), arrows = TRUE) #Plotting the inferred network
  
  setwd("C:/Users/rober/OneDrive/Desktop/Level - 3 networks Newick Notation/3a") #Where we are saving the file
  NN=read.evonet(folder[j])
  NN$edge.length = rep(1, dim(NN$edge)[1])
  plot(NN,arrows = T) #Plotting the true network for comparison
  
}










