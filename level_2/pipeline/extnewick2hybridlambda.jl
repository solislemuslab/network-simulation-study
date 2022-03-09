#GAB we need to modify this script so that it has two arguments: the file name of the extnewick which comes from SiPhyNetwork, and a filename for conversion to hybridlambda
#GAB Note that this script should be then operating over one single file at a time: One in, one out. The R script generate_datasets will take care of iterations
r_net_path = ARGS[1]

cd(r_net_path)#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/RNetworks")
using PhyloNetworks, PhyloPlots
using RCall

#net0 = readTopology("RNetwork")
#Nets=split(read(`ls`, String))
Nets=filter(x->startswith(x, "RNetwork"), readdir())

for i in 1:length(Nets)
	RNet = Nets[i]
	net0 = readTopology(RNet)
	#namefileJu = replace(RNet,"RNetwork"=>"RJuliaNet")*".txt"
	namefileHy = replace(RNet,"RNetwork"=>"JulHybrLamb")
	#cd(ju_net_path)#ju
	#writeTopology(net0,namefileJu)
	#cd(hyla_format_path)#hyla
	HyLaF = hybridlambdaformat(net0)
	open(namefileHy, "w") do io
		write(io, HyLaF)
	end
	#cd(r_net_path)
end 
