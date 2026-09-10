Goal: This simulation will asses the robustness of SNaQ when the true network is not level-1.

README: 

In order to run the simulation, you need the run following packages in R. (Note: some may already be pre-installed). This simulation will also require Julia.  

library(phangorn)
library(stringr)
library(MCMCprecision)
library(MSCquartets)



Lastly, each network was sketched using PhyloSketch, which can be found here: 

https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/phylosketch/

OpenSource Package: https://github.com/husonlab/phylosketch2  

Important: Each network from Phylosketch must be exported in Newick notation!(NEWICK file), which is a txt file after exporting. 



First, part of the simulation requires assigning a metric for all networks that are written in Newick notation imported in R. For each network constructed, they were categorized by their respective sample size, followed by their generators. Note, each file exported from PhyloSketch must be in a txt file. 
Call in a function that takes all networks (written in Newick Notation as a NEWICK file) and assign each edge a metric with gamma distribution (Shape and scale =1) and each hybridization parameter with a uniform distribution with values between (0.2,0.8). (Can be modified for high levels or based on preference.) 

Note that each function exported from Phylosketch will look something similar like this:
 ((((((((E)N#H3,F)Q)T#H2)O#H1,(((N#H3,D)L,C)K,B)J)I,A)M,(H,(O#H1,(T#H2,G)S)R)P)r;). 
  
											BELOW IS THE SIMULATION PIPELINE

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```r

#The goal is to remove the characters before the hybridization nodes and replace them with LGT followed by their respective parameters.

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
  
  #Assigning metric for each edge/ hybrid edge
  #Adding branch lengths
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
  } #adding :: and replacing Z for LGT# and the number of hybrid associated to it. Can be done for any level-k network.
    Y
}

```

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Can run this example in R to check to see if it works: 

#Example
X ="(((((((E)N#H3,F)Q)T#H2)O#H1,(((N#H3,D)L,C)K,B)J)I,A)M,(H,(O#H1,(T#H2,G)S)R)P)r;"
Y = top2metric(X)
Y


The result should return:

"(((((((E:0.642127725061995)#LGT3:2.50286903678723::0.724374739220366,F:3.50448497361176)Q:2.08678945904309)#LGT2:0.786253938458735::0.532265525078401)#LGT1:0.529362756536189::0.376071621477604,(((#LGT3:1.22191415703027::0.275625260779634,D:1.17560079074019)L:0.0794053253821385,C:1.69521327650384)K:0.109980137159687,B:2.1370628587849)J:0.308005926654248)I:0.757562118569836,A:0.195231384723383)M:0.0560168601245055,(H:1.49796833378233,(#LGT1:1.41698843570757::0.623928378522396,(#LGT2:0.136076517826142::0.467734474921599,G:0.0640945419308596)S:1.2778419290938)R:0.943117379529641)P:1.30563739388502)r;"



Next, for each network, we will simulate gene trees using PhyloCoalSimulator by the following pipeline in R: 
Note: The master file for simulating gene trees and inferring networks using SNaQ must have a the file changed from filename.R to filename.jl.
Below is the code for simulating gene trees 
										
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```r

#Just to clean the code up
path ='C:/Users/rober/OneDrive/Desktop/Level - 3 networks Newick Notation/'
path2master = 'C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/Level-3 Master File/masterfilelevel35000.jl'
path2masterSNaQ = 'C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/Level-3 Master File/SNaQMasterlevel35000.jl'


####Note this loop is for all networks in the semi-directed case
#######Sampling with the normal distribution from a folder network (with the correct path)

folder = list.files(paste0(path,"3a"))
for(j  in 1:length(folder)){
  file = folder[j]
  X = readLines(paste0(path, "3a/",file))
  X = paste0(X, collapse ="" )
  Y =removelab(X)
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


```
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



Next, inferring each sample of gene trees using SNaQ (note: hmax=1):


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```r

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

```
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Finally, we do a visual analysis between the inferred network and the true network below. 



--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```r

for(j  in 1:length(folder)){
  fileSNAQ =  paste0("C:/Users/rober/OneDrive/Desktop/R_Project_Folder_Phylogenetics/SNaQnetwork_Level3_5000/Generator3a/",folder[j],'snaqrun.networks')
  
  #####Goal of this is to read the text file, removes all the wording except to leave it in Newick Notation
  lineSNAQ = readLines(fileSNAQ)[1] #Reading the first line in the networks file
  lineSNAQ = gsub("\\, .*","",lineSNAQ) #Removes all wording except Newick notation
  lineSNAQ = gsub(":::",":1:",lineSNAQ) #Removes the triple colon and replaces it with only 1 colon followed by the branch length of 1
  networkSNAQ = read.evonet(text=lineSNAQ)
  networkSNAQ$edge.length =  rep(1, length(networkSNAQ$edge.length)) #Adding edge lengths to each remaining taxa, length 1
  plot(networkSNAQ, main = paste0("network SNaQ", j), arrows = TRUE)
  
  setwd("C:/Users/rober/OneDrive/Desktop/Level - 3 networks Newick Notation/3a")
  NN=read.evonet(folder[j])
  NN$edge.length = rep(1, dim(NN$edge)[1])
  plot(NN,arrows = T)
  
}

			
```
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

While the top is the simulation pipeline using R, this also requires the use of Julia for inferring each network and simulating gene trees: 

Below is simulating gene trees using Julia: 

Create a master file where R will call in the simulation pipeline. Below are the package required to run in Julia. 

										Below is the master file for both simulating gene trees and inferring them using SNaQ
													{filename.jl}
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```r

using PhyloNetworks
using CSV
using PhyloCoalSimulations
using DataFrames
using PhyNEST

cd("C:/Users/rober/OneDrive/Desktop/R Project Folder Phylogenetics/Thesis_Julia_Project_Level3_5000/3agenetrees") ### select where you want to save your files

num_gts = 5000; 
print("Number of gene trees is ", num_gts) 
graph_name = "testleveltwov4"  
gt_fname1 = string("gtrees_",graph_name) 
newick = "((((B:1,(C:1,(D:1,(E:1)#H2:1::0.5)L:1)J:1)H:1)#H1:1::0.5,A:1)G:1,(#H1:1::0.5,(#H2:1::0.5,F:1)O:1)M:1)K:1;"
ntwk = readTopology(newick);
gts1 = simulatecoalescent(ntwk, num_gts, 1; inheritancecorrelation=0.0 );
writeMultiTopology(gts1,gt_fname1)

```
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




										Below is inferring each network using SNaQ: 
												{filename.jl}


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```r

using PhyloNetworks;
gt_location = string("C:/Users/rober/OneDrive/Desktop/R Project Folder Phylogenetics/Thesis_Julia_Project_Level3_5000/",gt_fname1) 
cd("C:/Users/rober/OneDrive/Desktop/R Project Folder Phylogenetics/Thesis_Julia_Project_Level3_5000") #Where you want the file to be saved
using CSV;

trees = joinpath(gt_location); #From the previous gene tree you saved
genetrees = readMultiTopology(trees);
CF = readTrees2CF(trees);
q,t = countquartetsintrees(genetrees); 
df = writeTableCF(q,t)

CSV.write("tableCF.csv", df);
CF = readTableCF("tableCF.csv")
starttree = genetrees[2]
net1 = snaq!(starttree,CF, hmax=1, filename="net1", seed=1234)

```





