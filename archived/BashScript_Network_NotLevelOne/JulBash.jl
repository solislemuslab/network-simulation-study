#julia Runing bash
###
using PhyloNetworks, PhyloPlots
using RCall
using CSV

patFactors = "/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/CFactors/"
patGT = "/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/GeneTrees/"
patMaxCu = "/home/acosta/PhyloNetworks.jl.wiki/data_results/scripts"
patSnaq = "/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/SnaqNet/"

cd(patGT)
GenTrees = readdir()

#for i in 1:length(GenTrees)
for i in 1:3
	cd(patGT)
	GT1=GenTrees[i]
	TreSel = readMultiTopology(GT1)
	n1=replace(GT1, "_coal_unit" => ".csv")
	nameCF=replace(n1,"GeneTree"=>"CF")
	n2=replace(GT1, "_coal_unit" => "")
	nameSnaq=replace(n2,"GeneTree"=>"SnaqNet")
	cd(patFactors)
	q,t = countquartetsintrees(TreSel;)
	df = writeTableCF(q,t)
	CSV.write(nameCF, df);
	cd(patMaxCu)
	run(Cmd(["../scripts/get-pop-tree.pl",patFactors*nameCF]))
	cd(patFactors)
	raxmlCF = readTableCF(nameCF)
	namTree=nameCF*".QMC.tre"
	TreeSnaq = readTopology(namTree)
	#
	net0 = snaq!(TreeSnaq,  raxmlCF, hmax=0, seed=123,runs=2)
	cd(patSnaq)
	writeTopology(net0,nameSnaq)
end
