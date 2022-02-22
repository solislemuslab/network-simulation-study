##
using PhyloNetworks, PhyloPlots, RCall

UltRNetwork = "((((t6:0.7351678649,t8:0.7351678649):0.5285122109,((t2:0.3278540562,#H18:0::0.7669243112):0.9358260196)#H13:0::0.8162281318):0.9866902289,(((t4:0.05074118239,t3:0.05074118239):0.6202557755,(((t5:0.1598060242,t1:0.1598060242):0.1047496504,t7:0.2645556746):0.06329838161)#H18:0.3431429017::0.2330756888):0.5926831179,#H13:0::0.1837718682):0.9866902289):1.31728298);"


net0 = readTopology(UltRNetwork)
plot(net0, :R, showEdgeNumber=true,showNodeNumber=true,showTipLabel=true);
####


function NetChar(net0)
	Edge1 = net0.edge
	Len = Vector{Float64}()
	IsLeaf = Vector{Bool}()
	IsHyb = Vector{Bool}()
	Parent = Vector{Int}()
	Child = Vector{Int}()
	Id = Vector{Int}()
	Tips = Vector{Int}()
	#IsMaj = Vector{Bool}()

	for i in 1:length(Edge1)
		ed1 = Edge1[i]
		ismaj = ed1.isMajor

		len1 = ed1.length
		#hyb = ed1.hybrid
		ch1 = ed1.node[1]
		lef1 = ch1.leaf
		childN = ch1.number
		parentN = ed1.node[2].number

		if lef1
			push!(Tips,childN)
		end	

		#push!(Id,i)
		push!(Len,len1 )
		#push!(IsLeaf,lef1)
		push!(Parent,parentN)
		push!(Child,childN)
			#push!(IsHyb,hyb)

	end
	return (Len=Len,Parent=Parent,Child=Child,Tips=Tips)
end

function LengNodeRoot(Len,Parent,Child,Tip)
	ini=Tip
	indexCh1 = findall( x -> x == ini[1], Child)
	ini = Parent[indexCh1]
	L1 = Len[indexCh1]
	logi1 = length(indexCh1)==1

	LenTot = Vector{Float64}()
	push!(LenTot,L1[1])

	while logi1
		indexCh1 = findall( x -> x == ini[1], Child)
		if length(indexCh1)==2
			indexCh1 = findall( x -> x == ini[1], Child)[1]
		end	
		ini = Parent[indexCh1]
		L1 = Len[indexCh1]
		logi1 = length(indexCh1)==1
		if logi1
			push!(LenTot,L1[1])
		end
	end
	return sum(LenTot)
end



NetDes = NetChar(net0);


LengNodeRoot(NetDes.Len, NetDes.Parent, NetDes.Child, NetDes.Tips[3])





tp = NetDes.Tips

LTipRoot = Vector{Float64}()
for i in 1:length(tp)
	res = LengNodeRoot(NetDes.Len, NetDes.Parent, NetDes.Child, tp[i])
	push!(LTipRoot,res)
end

NotUlt=LTipRoot	














a=writeTopology(net0)
Net=a


indexNot0 = findall( x -> x >0, NetDes.Len)
BraLe = NetDes.Len[indexNot0]
MinEd = min(BraLe...)

Factor=(1/MinEd)*0.5

function scaleBranches!(Net,Factor)
    Net=a
	x1 = Vector{Int}()
	x2 = Vector{Int}()
	push!(x2,1)

	for i in 1:length(Net)-1
		if (Net[[i]]==":")&&(Net[[i+1]]!=":")
			push!(x1,i)
		end
		if (Net[[i]]==",")|| (Net[i:i+1]=="::")||(Net[[i]]==")")
			push!(x2,i)
		end	

	end

	push!(x1,length(Net))


	NumLen = Vector{Float64}()
	for j in 1:length(x1)-1
		#j=1
		num = parse(Float64, Net[x1[j]+1:x2[j+1]-1])
		push!(NumLen,num)
	end
	
	NumLenNew = NumLen*Factor

  
	NetRes=""
	for i in 1:length(NumLenNew)
		ro=Net[x2[i]:x1[i]]*string(NumLenNew[i])
		NetRes=NetRes*ro
	end
	NetResF = NetRes*");"

	return NetResF
end

rol=scaleBranches!(a,100)



net1 = readTopology(rol)



plot(net0, :R, showEdgeNumber=true,showNodeNumber=true,showTipLabel=true);

###

rol=scaleBranches!(a,100000)
net1 = readTopology(rol)
HyLaF = hybridlambdaformat(net1)
open("HL_3", "w") do io
	write(io, HyLaF)
end









