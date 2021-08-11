#!/usr/bin/julia
#print("hi Julia21rg")
#print("ayuda")
#
IniPat = pwd()
cd(joinpath((IniPat),"RNetworks"))
#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/RNetworks")
using PhyloNetworks, PhyloPlots
using RCall

<<<<<<< HEAD
Nets=split(read(`ls`, String))

=======
>>>>>>> e5fc8db3fd3d2181541ea303d6a023ce50021dde
#net0 = readTopology("RNetwork")
#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne")
#writeTopology(net0, "JuliaNetwork1.txt")
#R"pdf(file='NetworkjuliaR2.pdf')"
#plot(net0, :R)
#R"dev.off()"

nets_not_level1 = readMultiTopology("nets_not_level1.extnex")

# bulk format conversion of networks read in nets_not_level1
open("nets_not_level1.hybridlambda", "w") do io    
    for i = 1:length(nets_not_level1)
        hyblambFormatted = hybridlambdaformat(nets_not_level1[i])
        write(io, "$hyblambFormatted\n")
        flush(io)
    end
end
    

#savefig("plotJu.png")
#print(hybridlambdaformat(net0))
#cd("/home/acosta/Escritorio/Working/Filogenetica/Bash")
#readdir()
<<<<<<< HEAD

#



for i in 1:length(Nets)
	cd(joinpath((IniPat),"RNetworks"))
	RNet = Nets[i]
	net0 = readTopology(RNet)
	namefileJu = replace(RNet,"RNetwork"=>"RJuliaNet")*".txt"
	namefileHy = replace(RNet,"RNetwork"=>"JulHybrLamb")
	cd(joinpath((IniPat),"RJuliaNetworks"))
	writeTopology(net0,namefileJu)
	cd(joinpath((IniPat),"JuliaForHybLam"))
	HyLaF = hybridlambdaformat(net0)
	open(namefileHy, "w") do io
		write(io, HyLaF)
	end
	cd(joinpath((IniPat),"RNetworks"))
end
=======
>>>>>>> e5fc8db3fd3d2181541ea303d6a023ce50021dde
