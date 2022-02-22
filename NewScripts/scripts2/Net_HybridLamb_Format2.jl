

r_net_path = ARGS[1]
ju_net_path = ARGS[2]
hyla_format_path = ARGS[3]

cd(r_net_path)#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/RNetworks")
using PhyloNetworks, PhyloPlots
using RCall

#net0 = readTopology("RNetwork")
Nets=split(read(`ls`, String))


for i in 1:length(Nets)
	cd(r_net_path)
	RNet = Nets[i]
	net0 = readTopology(RNet)
	namefileJu = replace(RNet,"RNetwork"=>"RJuliaNet")*".txt"
	namefileHy = replace(RNet,"RNetwork"=>"JulHybrLamb")
	cd(ju_net_path)#ju
	writeTopology(net0,namefileJu)
	cd(hyla_format_path)#hyla
	HyLaF = hybridlambdaformat(net0)
	open(namefileHy, "w") do io
		write(io, HyLaF)
	end
	cd(r_net_path)
end 
