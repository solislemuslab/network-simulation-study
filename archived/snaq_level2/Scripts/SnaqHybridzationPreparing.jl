
#LIBRARIES
using PhyloNetworks, PhyloPlots
using RCall
using CSV
###
using DataFrames
using DelimitedFiles


# WHAT DO?
#get concordance factors
#get Conseonsous tree using qmct and concordace factors
#Storage it in suiatble path ready to run with snaq
##For each network a folder is created as "snaq_input_[i]" where [i] is the number of the network
## in there will be storage the Concordance factors, Consensous tree (needed for snaq) and other supplementary files



#ARGUMENTS
#For do that the following arguments are needed:
#This argumenst are obtained from "Species_Gene_Networks.sh"

#path of the concordance factors ("$PATH_ULT_C_FACTOR")
path_cfactors = ARGS[1]*"/"

#path of gene trees obtained with hybrid lambda ("$PATH_ULT_GT")
path_gtrees = ARGS[2]*"/"

#path of scripts of quarter maxcut ("$PTH_QMCT")
path_qmct = ARGS[3]*"/"

#pat of the trees obtained from quarter maxcut("$PTH_QMCT" )
path_tree_qmct = ARGS[4]*"/"

#path of networks obtained with R ("$PATH_ULT_R" )
path_rnetwork = ARGS[5]*"/"

#path of network in hybrid lambda format obtained with julia("$PATH_ULT_JHL")
ParHL = ARGS[6]*"/"

#Path to store the results as input for run snaq in the next step ("$PTH_SNAQ_INPUT")#path_store_snaq_input
path_store_snaq_input = ARGS[7]*"/"

#Path where there are stored the arguments used to create the networks an gene trees ("$PATH_ULT_README")
#This will considerer as supplementary files
path_readme = ARGS[8]*"/"




#START

#Get concordance factors and consesous tree with qmct
cd(path_gtrees)
name_genetrees = readdir()

for i in 1:length(name_genetrees)
	cd(path_gtrees)
	nam_gt1=name_genetrees[i]
	
	#Code to delet the sub indice "_1" that Hybrid Lambda add
	open_file = open(nam_gt1)
	file_lines = readlines(open_file)

	for line in 1:length(file_lines)
		file_lines[line] = replace(file_lines[line ], r"_1:" => s":")
	end
	###finish code that delet subindice	"_1"
	
	#get suitable name for Cconcordance factors using GT names
	nam_cf = replace(nam_gt1, "_coal_unit" => ".csv")
	name_cfactors = replace(nam_cf,"GeneTree"=>"CF")
	
	#Get the concordance factors
	gene_tree = readMultiTopology(nam_gt1)
	cd(path_cfactors)
	q,t = countquartetsintrees(gene_tree;)
	df = writeTableCF(q,t)

	
#Gustavo code to fix the pop tree omit	
	cfs = readTableCF(df)
	zeros = repeat([0], nrow(df))

	df_zeros = DataFrame(CF12_34_lo=zeros,
						 CF12_34_hi=zeros,
						 CF13_24_lo=zeros,
						 CF13_24_hi=zeros,
						 CF14_23_lo=zeros,
						 CF14_23_hi=zeros)

	df_ticr = hcat(select(df, 1:5),
				   select(df_zeros, 1:2),
				   select(df, 6),
				   select(df_zeros, 3:4),
				   select(df, 7),
				   select(df_zeros, 5:6))	
#Gustavo code end
	
	
	#save the concordance factors
	CSV.write(name_cfactors, df_ticr)
	
	#Run cuarter max cut
	cd(path_qmct)
	run(Cmd(["../scripts/get-pop-tree.pl",path_cfactors*name_cfactors]))
	
	cd(path_cfactors)
	pat1 = path_cfactors*name_cfactors*".QMC.tre"
	pat2 = path_tree_qmct*name_cfactors*".QMC.tre"
	mv(pat1,pat2)
end



#Transfer the "conseosous tree", "concordance factors" and all the supplementary files to the path snaq_input
#we use the name of r networks to nominate the folder as "snaq_input_[i]"
cd(path_rnetwork)
name_rnetworks = readdir()

for i in 1:length(name_rnetworks)
	cd(path_store_snaq_input)
	NumPath = replace(name_rnetworks[i], "RNetwork_" => "")
	mkdir("snaq_input_"*NumPath)
	
	#transfer suplementary files to the path snaq_input
	PatOriR = path_rnetwork*"RNetwork_"*NumPath
	PatDesR = path_store_snaq_input*"snaq_input_"*NumPath*"/"*"RNetwork_"*NumPath
	cp(PatOriR,PatDesR)
	
	PatOriHL = ParHL*"JulHybrLamb_"*NumPath
	PatDesHL = path_store_snaq_input*"snaq_input_"*NumPath*"/"*"JulHybrLamb_"*NumPath
	cp(PatOriHL,PatDesHL)
	
	PatOriGT = path_gtrees*"GeneTree_"*NumPath*"_coal_unit"
	PatDesGT = path_store_snaq_input*"snaq_input_"*NumPath*"/"*"GeneTree_"*NumPath*"_coal_unit"
	cp(PatOriGT,PatDesGT)
	
	#transfer main files (CF and consensous tree) to the path snaq_input
	
	PatOriCF = path_cfactors*"CF_"*NumPath*".csv"
	PatDesCF = path_store_snaq_input*"snaq_input_"*NumPath*"/"*"CF_"*NumPath*".csv"
	cp(PatOriCF,PatDesCF)
	
	PatOriCMQ = path_tree_qmct*"CF_"*NumPath*".csv.QMC.tre"
	PatDesCMQ = path_store_snaq_input*"snaq_input_"*NumPath*"/"*"CF_"*NumPath*".csv.QMC.tre"
	cp(PatOriCMQ,PatDesCMQ)
	
	PatOriReadme= path_readme*"readme_"*NumPath*".txt"
	PatDesReadme = path_store_snaq_input*"snaq_input_"*NumPath*"/"*"readme_"*NumPath*".txt"
	cp(PatOriReadme,PatDesReadme)
	
end

