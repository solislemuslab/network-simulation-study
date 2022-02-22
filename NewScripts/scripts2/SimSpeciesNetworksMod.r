
# 1. Read all the parameters needed to do the simmulation
# 2. Create all the paths needed
# 	 2.1. Output path to save the main results (Good RUltrametric Networks and his Gene trees)
#	 2.2. Temporary path to save the first results to test "complexity" (Extrange Networks) and ultrametricity
# 3. Get Networks without warnings: networks (avoid trees) without warnings
#	 3.1. Simulate networks
#	 3.2. Ommit Extrange Networks 
#	 3.3. Ommit Trees (Only select networks)
#	 3.4. Save all the good Networks
# 4. Select Ultrametric networks
# 	 4.1. Convert the good networks in hybrid lambda format using julia script ("Net_HybridLamb_Format2.jl")
#    4.2. Copy the good networks in hybrid lambda path, run it and save the ultrametricity results.
#    4.3. Read The ultrametricity results to know if is or not Ultrametric and save where it belongs
# 5. Get Gene Trees using Hybrid lambda 

#####################################################################################
#1. Parsmeters of the simulation
##################################################################3
#numbsim1=100
#ntips=8
#lambda=0.9
#mu=0
#nu=0.09
#hybprops=c(0.5,0.25,0.25)
#num_gene_tree<-500
#genseed<-25

numbsim1 <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
ntips <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
lambda <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
mu <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
nu <- as.numeric(commandArgs(trailingOnly=TRUE)[5])
hybpropsr <- commandArgs(trailingOnly=TRUE)[6]
num_gene_tree <- as.numeric(commandArgs(trailingOnly=TRUE)[7])
genseed <- as.numeric(commandArgs(trailingOnly=TRUE)[8])

######################################################################################
eval(parse(text=hybpropsr))
hybprob <- paste0(hybprops,collapse="-")
hyblam_path <-"/home/carlos/hybrid-Lambda-0.6.2-beta/src"#Hibrid lambda src path (fixed), it could be considerer as parameter






#####################################################
##2. Creating all the paths
#############################################
scrip_path <- getwd()#path where are all the scripts needed
setwd('..')#Go to the parent path
parent_pat <- getwd()#Parent path

#Temporaty path 
temp_path <- paste(parent_pat,"/","temporary_path",sep="")

#Into the temporary path save rnetwork, rextrange network, 
temp_r_net <- paste(temp_path ,"/","r_network",sep="")
temp_rextrange_net <- paste(temp_path ,"/","r_extrange_network",sep="")
temp_ju_net <- paste(temp_path ,"/","jul_network",sep="")
temp_hybla_net <- paste(temp_path ,"/","hyla_network",sep="")

# Also in temporary path save the ultramecric network
temp_ult_met <- paste(temp_path ,"/","Ultrametric",sep="")
temp_ult_met_hylaformat <- paste(temp_ult_met,"/","hyla_format",sep="")

temp_not_ult_met <- paste(temp_path ,"/","NotUltrametric",sep="")
temp_not_ult_met_hylaformat <- paste(temp_not_ult_met,"/","hyla_format",sep="")
#In output pat save only the genetrees.
out_path<- paste(parent_pat,"/","Output",sep="")

dir.create(path=temp_path)
dir.create(path=temp_r_net)
dir.create(path=temp_rextrange_net)
dir.create(path=temp_hybla_net)
dir.create(path=temp_ju_net)
dir.create(path=temp_ult_met)
dir.create(path=temp_not_ult_met)
dir.create(path=out_path)
dir.create(path=temp_ult_met_hylaformat)
dir.create(path=temp_not_ult_met_hylaformat)
######################################################################3





#########################################
# Function to identify warnings and use that to identify the extrange networks
###############################################
#the following function help me to manage the warnings and I found it in this link:
#https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 
###############################################################################33





###################################################
#3. Get Networks without warnings
###################################################
library(SiPhyNetwork)#library to simulate Networks, this use "ape" as dependence
library(geiger)#is.phylo()##Information about tree goes extinct=0 and no extinct tips are sampled=1

#3.1 Simulate networks
set.seed(genseed)
bdNS1 <- sim.bdh.taxa.ssa(n=ntips, 
                                 numbsim = numbsim1, 
                                 lambda = lambda, 
                                 mu = mu , 
                                 nu = nu,
                                 hybprops = hybprops, 
                                 hyb.inher.fxn = make.beta.draw(1, 1), 
                                 frac = 1,
                                 mrca = FALSE, 
                                 complete = TRUE, 
                                 stochsampling = FALSE, 
                                 hyb.rate.fxn = NULL, 
                                 trait.model = NULL)
#
bdNS1 <- bdNS1[!sapply(X = bdNS1, FUN = is.null)]
#Information about tree goes extinct=0 and no extinct tips are sampled=1
bdNS1 <- bdNS1[sapply(X = bdNS1, FUN = is.phylo)]





##############################
#3.2. Ommit Extrange Networks
##############################
setwd(temp_r_net)#path
outWarn<-c()
for(i in 1:length(bdNS1)){
tree<-bdNS1[[i]]
warnn1<-length(catchToList(plot(tree,main="R Network",cex=1))$warnings)
outWarn<-c(outWarn,warnn1)
}
bdNS<-bdNS1[which(outWarn==0)]#Good Networks 

#The procces produce pdf, delete it because it generate fatal error with Julia.
x<-dir()
file.remove(x)



##########################################
#3.3. Ommit Trees (Only select networks)
#######################################
out<-c()
for(i in 1:length(bdNS)){
	tre_i <- bdNS[[i]]
	ret_numb <- nrow(tre_i$reticulation)
	res <-ret_numb>=1
	out<-c(out,res)
}


##################################3
#3.4. Save all the good Networks
####################################3
id<-which(out)
for(j in 1:length(id)){
tree<-bdNS[[id[j]]]
write.net(tree,file=paste("RNetwork_",j,sep=""))
}


#####################################
#Save the extrange networks
#####################################
setwd(temp_rextrange_net)#Goto path of extrange networks
bdNSWarnings<-bdNS1[which(outWarn>0)]
for(k in 1:length(bdNSWarnings)){
tree<-bdNSWarnings[[k]]
write.net(tree,file=paste("RNetworkWarn_",k,sep=""))
}




################################################################################
# 4. Select Ultrametric networks
#####################################################################

#4.1 convert the good networks in hybrid lambda format using julia
setwd(scrip_path) #Go to the path of scripts and run julia code to 
ju_script <- paste("julia Net_HybridLamb_Format2.jl", temp_r_net, temp_ju_net, temp_hybla_net)
system(ju_script)


##############################################################
#4.2 Copy networks in Hybrid lambda format to Hybrid lambda/src to run Hybrid lambda
list_of_files <- list.files(temp_hybla_net) 
file.copy(file.path(temp_hybla_net,list_of_files), hyblam_path)



################################################################
#4.3. Read The ultrametricity results to know if is or not Ultrametric and save where it belongs
setwd(hyblam_path)# go to Hybrid lambda/src path

for(i in 1:length(list_of_files)){

hyla_network_i <- list_of_files[i]
#Extrac subindice to rename the following outputs
sub_i <- substr(x =hyla_network_i, start=12, stop=nchar(hyla_network_i))

NameforHyLambda<-paste("HybUlt",sub_i,sep="")
name_out <- paste("outerr",sub_i,sep="")
	
#Run hybrid lambda using R
hylamb_instructions <- paste("../src/hybrid-Lambda -spcu ../src/", hyla_network_i, " -dot -label -o ", NameforHyLambda, " > ", name_out, ".txt 2>&1",sep="")
system(hylamb_instructions)

filename_ultrametricity <- paste(name_out, ".txt",sep="")
	
ultrametric_logic <- grep(x=readLines(filename_ultrametricity)[3],pattern="Dot figure generated in file")
if(length(ultrametric_logic)==0){
	#Not Ultrametric, we can add more ultrametric outputs
	file.copy(paste(hyblam_path,"/",hyla_network_i,sep=""), temp_not_ult_met_hylaformat)#hybrid lambda format
}
	else {
		#Ultrametric, we can add more Not ultrametric outputs
		file.copy(paste(hyblam_path,"/",hyla_network_i,sep=""), temp_ult_met_hylaformat)#hybrid lambda format
}
}


#Delete all files added in Hybrid lambda/src
all_files_add_hyla<-c(dir(pattern="HybUlt"), dir(pattern="outerr"), dir(pattern="JulHybrLamb_"))
file.remove(file.path(hyblam_path,all_files_add_hyla))








#########################################
#5. Get Gene Trees using Hybrid lambda 
#########################################
#(../src/hybrid-Lambda -spcu ../src/JuliaHibrid2 -num 500 -seed 2 -o example3)
#Copy all Ultrametric networks to hybrid lambda to get gene trees

ultrametric_networks <- dir(path=temp_ult_met_hylaformat)
file.copy(paste(temp_ult_met_hylaformat,"/",ultrametric_networks,sep=""), hyblam_path)


for( i in 1:length(ultrametric_networks)){
	ult_net_i <- ultrametric_networks[i]
	
	Net<-gene_tree_name <- gsub(pattern="JulHybrLamb", replacement="Net", x=ultrametric_networks[i])
	hyla_seed <-sample(1:10000,1)

	nam_path_out_Ni <- paste(Net,"_numbsim1-", numbsim1, "_ntips-", ntips,"_lambda-" ,lambda,"_mu-", mu,"_nu-", nu,"_hybprops-",hybprob ,
							 "_GTnum-",num_gene_tree,"_Rseed_",genseed,"_hylaseed_",hyla_seed,sep="")
		path_out_i <- paste(out_path,"/",nam_path_out_Ni,sep="")
	path_out_i 
	dir.create(path=path_out_i)
	
	gene_tree_name <- gsub(pattern="JulHybrLamb", replacement="GT", x=ultrametric_networks[i])
	hylamb_instructions_genetree <- paste("../src/hybrid-Lambda -spcu ../src/", ult_net_i, " -num ", num_gene_tree, " -seed ", hyla_seed ," -o ", gene_tree_name,sep="")
	hylamb_instructions_genetree
	system(hylamb_instructions_genetree)
	#Copy the gene trees in the output path
	file.copy(paste(hyblam_path,"/",paste(gene_tree_name,"_coal_unit",sep=""),sep=""), path_out_i)
	
	#Copy RUltrametric Networks in the output path
	rnet_name <- gsub(pattern="JulHybrLamb", replacement="RNetwork", x=ultrametric_networks[i])
	file.copy(paste(temp_r_net,"/", rnet_name,sep=""), path_out_i)
}



#Delete all files added in Hybrid lambda/src
all_files_add_hyla<-c(dir(pattern="GT"), dir(pattern="JulHybrLamb_"))
file.remove(file.path(hyblam_path,all_files_add_hyla))



