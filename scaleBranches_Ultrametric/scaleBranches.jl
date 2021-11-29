##
using PhyloNetworks, PhyloPlots, RCall



UltRNetwork = "((((t6:0.7351678649,t8:0.7351678649):0.5285122109,((t2:0.3278540562,#H18:0::0.7669243112):0.9358260196)#H13:0::0.8162281318):0.9866902289,(((t4:0.05074118239,t3:0.05074118239):0.6202557755,(((t5:0.1598060242,t1:0.1598060242):0.1047496504,t7:0.2645556746):0.06329838161)#H18:0.3431429017::0.2330756888):0.5926831179,#H13:0::0.1837718682):0.9866902289):1.31728298);"

net0 = readTopology(UltRNetwork)


plot(net0, :R, showEdgeNumber=true,showNodeNumber=true,showTipLabel=true);


NetJuliaFromat = writeTopology(net0)


#Function that scale the branches lengths multiplying by a factor
function scaleBranches!(Net,Factor)
    x1 = Vector{Int}()#Vetor to store the position when start one branch length (after ":")
	x2 = Vector{Int}()#Vetor to store the position when finish one branch length (before "::" or ")")
	push!(x2,1)

	for i in 1:length(Net)-1
		if (Net[[i]]==":")&&(Net[[i+1]]!=":")# after ":" is posible to find a number or one more time ":", we need to find a number
			push!(x1,i)
		end
		if (Net[[i]]==",")|| (Net[i:i+1]=="::")||(Net[[i]]==")")
			push!(x2,i)
		end	
	end

	push!(x1,length(Net))

	#Extract all the banches lengths 
	NumLen = Vector{Float64}()
	for j in 1:length(x1)-1
		num = parse(Float64, Net[x1[j]+1:x2[j+1]-1])
		push!(NumLen,num)
	end
	
	NumLenNew = NumLen*Factor#Multiply by a factor
	
	#Reconstructing the network with the new lengths
	NetRes=""
	for i in 1:length(NumLenNew)
		ro=Net[x2[i]:x1[i]]*string(NumLenNew[i])
		NetRes=NetRes*ro
	end
	NetResF = NetRes*");"

	return NetResF
end


NetScaled = scaleBranches!(NetJuliaFromat,1000)
NetScaled





NetScaledJ = readTopology(NetScaled)
NetHyLaF = hybridlambdaformat(NetScaledJ)
open("HL_3", "w") do io
	write(io, NetHyLaF)
end



###########################################################################################################################3
##cd /home/carlos/hybrid-Lambda-0.6.2-beta/src
#../src/hybrid-Lambda -spcu ../src/NetworkNameForHLFormat -dot -label -o NameforHyLambda > FilenameToStoreOutput.txt 2>&1
##../src/hybrid-Lambda -spcu ../src/HL_3 -dot -label -o Net3 > FilenameToStoreOutput.txt 2>&1
