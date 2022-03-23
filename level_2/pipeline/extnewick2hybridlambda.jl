rnetwork=ARGS[1]
#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/RNetworks")
using PhyloNetworks

a=split(rnetwork,"_")

netini = readTopology(a[1])
HyLaFini = hybridlambdaformat(netini)
print(HyLaFini)

for i in 2:length(a)
	print("_")
	net0 = readTopology(a[i])
	HyLaF1 = hybridlambdaformat(net0)
	#HyLaF0 = HyLaF0*"_"*HyLaF1
	print(HyLaF1)
	
end
