
###
using PhyloNetworks, PhyloPlots
using RCall
using CSV
using DataFrames

#patFactors = "/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/BashSnaq/Output/SnaqOut/CFactors/"
#patFactors = "/home/carlos/Testing/Filogenetica/Test1/Output/UltraMetric/CFactors/"
patCFactors = ARGS[1]*"/"
#patGT = "/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/BashSnaq/Output/SnaqOut/GeneTrees/"
#patGT = "/home/carlos/Testing/Filogenetica/Test1/Output/UltraMetric/GeneTrees/"
patGT = ARGS[2]*"/"
#patMaxCu = "/home/acosta/PhyloNetworks.jl.wiki/data_results/scripts"#FixedPath
#patMaxCu = "/home/carlos/phyloNet/PhyloNetworks.jl.wiki/data_results/scripts/"#FixedPath
patMaxCu = ARGS[3]*"/"
#patSnaq = "/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/BashSnaq/Output/SnaqOut/SnaqNet/"

#pat to storage the tree from maxcut
#TreMaxCu = "/home/carlos/Testing/Filogenetica/Test1/Output/UltraMetric/QMCtree/"
TreMaxCu = ARGS[4]*"/"


PatRnet = ARGS[5]*"/"
ParHL = ARGS[6]*"/"

patStorSnaq = ARGS[7]*"/"

PatReadme = ARGS[8]*"/"

cd(patGT)
GenTrees = readdir()

for i in 1:length(GenTrees)
	cd(patGT)
	GT1=GenTrees[i]
	TreSel = readMultiTopology(GT1)
	n1=replace(GT1, "_coal_unit" => ".csv")
	nameCF=replace(n1,"GeneTree"=>"CF")
	n2=replace(GT1, "_coal_unit" => "")
	nameSnaq=replace(n2,"GeneTree"=>"SnaqNet")
	cd(patCFactors)
	q,t = countquartetsintrees(TreSel;)
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
	
	CSV.write(nameCF, df_ticr)
	
	cd(patMaxCu)
	run(Cmd(["../scripts/get-pop-tree.pl",patCFactors*nameCF]))
	
	cd(patCFactors)
	Pat1 = patCFactors*nameCF*".QMC.tre"
	Pat2 = TreMaxCu*nameCF*".QMC.tre"
	mv(Pat1,Pat2)
end


cd(PatRnet)
NetSnaq = readdir()

for i in 1:length(NetSnaq)
	cd(patStorSnaq)
	NumPath = replace(NetSnaq[i], "RNetwork_" => "")
	mkdir("SnPatt_"*NumPath)
	
	PatOriR = PatRnet*"RNetwork_"*NumPath
	PatDesR = patStorSnaq*"SnPatt_"*NumPath*"/"*"RNetwork_"*NumPath
	cp(PatOriR,PatDesR)
	
	PatOriHL = ParHL*"JulHybrLamb_"*NumPath
	PatDesHL = patStorSnaq*"SnPatt_"*NumPath*"/"*"JulHybrLamb_"*NumPath
	cp(PatOriHL,PatDesHL)
	
	PatOriGT = patGT*"GeneTree_"*NumPath*"_coal_unit"
	PatDesGT = patStorSnaq*"SnPatt_"*NumPath*"/"*"GeneTree_"*NumPath*"_coal_unit"
	cp(PatOriGT,PatDesGT)
	
	PatOriCF = patCFactors*"CF_"*NumPath*".csv"
	PatDesCF = patStorSnaq*"SnPatt_"*NumPath*"/"*"CF_"*NumPath*".csv"
	cp(PatOriCF,PatDesCF)
	
	PatOriCMQ = TreMaxCu*"CF_"*NumPath*".csv.QMC.tre"
	PatDesCMQ = patStorSnaq*"SnPatt_"*NumPath*"/"*"CF_"*NumPath*".csv.QMC.tre"
	cp(PatOriCMQ,PatDesCMQ)
	
	PatOriReadme= PatReadme*"readme_"*NumPath*".txt"
	PatDesReadme = patStorSnaq*"SnPatt_"*NumPath*"/"*"readme_"*NumPath*".txt"
	cp(PatOriReadme,PatDesReadme)
	
end


