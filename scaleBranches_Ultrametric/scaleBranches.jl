using PhyloNetworks, PhyloPlots, RCall

"""
Code
"""

#Function that scales branch lengths multiplying by a factor

#GAB avoid capitalised variables, these are best restricted to classes rather than variables.

#GAB The ! functions modify their input, therefore it is unnecessary to assign them to another object. Both netJuliaFormat and netScaled would have the same content. The alternative is to avoid the ! syntax in this function. In fact, the function below runs perfectly well without the ! operator

#GAB document your arguments:
# net, a network of class XXXXX
# factor, an integer for which all the edges will be multiplied

function scaleBranches(net, factor)
    x1 = Vector{Int}() # Vector to store the position when start one branch length (after ":")
    x2 = Vector{Int}() # Vector to store the position when finish one branch length (before "::" or ")")
    push!(x2,1)

    # explain in words what your code does so that the reviewer can assess whether it just does what you intended to
    for i in 1:length(net)-1
        if (net[[i]]==":")&(net[[i+1]]!=":")# after ":" is posible to find a number or one more time ":", we need to find a number
	    push!(x1,i)
	end
	if (net[[i]]==",")|(net[i:i+1]=="::")|(net[[i]]==")")
	    push!(x2,i)
	end	
    end

    push!(x1,length(net))

    #Extract all the branch lengths 
    numLen = Vector{Float64}()
    for j in 1:length(x1)-1
	num = parse(Float64, net[x1[j]+1:x2[j+1]-1])
	push!(numLen,num)
    end
	
    numLenNew = numLen*factor#Multiply by a factor
	
    #Reconstructing the network with the new lengths
    netRes=""
    for i in 1:length(numLenNew)
	ro=net[x2[i]:x1[i]]*string(numLenNew[i])
	netRes=netRes*ro
    end
    netResF = netRes*");"

    return netResF
end

"""
Testing the function
"""

# Try to use first small examples for testing your code. This net is both too complex and too large to spot coding errors in case they happen. I show below a tiny example on which your function works, and looks fine.
UltRNetwork = "((((t6:0.7351678649,t8:0.7351678649):0.5285122109,((t2:0.3278540562,#H18:0::0.7669243112):0.9358260196)#H13:0::0.8162281318):0.9866902289,(((t4:0.05074118239,t3:0.05074118239):0.6202557755,(((t5:0.1598060242,t1:0.1598060242):0.1047496504,t7:0.2645556746):0.06329838161)#H18:0.3431429017::0.2330756888):0.5926831179,#H13:0::0.1837718682):0.9866902289):1.31728298);"

net0 = readTopology(UltRNetwork)

plot(net0, :R, showEdgeNumber=true,showNodeNumber=true,showTipLabel=true)

#GAB avoid capitalised variables, these are best restricted to classes rather than variables.
0netJuliaFormat = writeTopology(net0)

#GAB The code below tests the function
netScaled = scaleBranches!(netJuliaFormat,1000)
netScaled = readTopology(netScaled) # we need a network object, not a string
netScaled

#GAB the plot below looks adequate, at least the skeleton looks the same
plot(netScaled, :R, showEdgeNumber=true,showNodeNumber=true,showTipLabel=true)

""" 
Unit testing
"""

#GAB Test 1: A small tree with two leaves and unity branch lengths

miniNet = "(t1:1,t2:1);"
#miniNet = readTopology(miniNet)
plot(miniNet, :R, showEdgeNumber=true,showNodeNumber=true,showTipLabel=true)
miniNet = scaleBranches(miniNet,1000)
#GAB Are branch lengths equal to 1000? Yes, they are, congrats!


#GAB Test2: Multiplying by unity leaves the branch lengths unchanged
#GAB Good job! your function leaves unchanged the testing small tree!
miniNet = "(t1:1.0,t2:1.0);"
miniNet == scaleBranches(miniNet, 1) # this is true
miniNet === scaleBranches(miniNet, 1) # this is true

#GAB Bad news for the largest network :(
net0 = "((((t6:0.7351678649,t8:0.7351678649):0.5285122109,((t2:0.3278540562,#H18:0::0.7669243112):0.9358260196)#H13:0::0.8162281318):0.9866902289,(((t4:0.05074118239,t3:0.05074118239):0.6202557755,(((t5:0.1598060242,t1:0.1598060242):0.1047496504,t7:0.2645556746):0.06329838161)#H18:0.3431429017::0.2330756888):0.5926831179,#H13:0::0.1837718682):0.9866902289):1.31728298);"
net00 = "((((t6:0.7351678649,t8:0.7351678649):0.5285122109,((t2:0.3278540562,#H18:0.0::0.7669243112):0.9358260196)#H13:0.0::0.8162281318):0.9866902289,(((t4:0.05074118239,t3:0.05074118239):0.6202557755,(((t5:0.1598060242,t1:0.1598060242):0.1047496504,t7:0.2645556746):0.06329838161)#H18:0.3431429017::0.2330756888):0.5926831179,#H13:0.0::0.1837718682):0.9866902289):1.31728298);"
#GAB the two networks should be identical because they are multiplied by 1.0; however, they are not.
net0 == scaleBranches(net0, 1.0) # this is false
net0 === scaleBranches(net0, 1.0) # this is false
net00 == scaleBranches(net00, 1.0) # this is true
net00 === scaleBranches(net00, 1.0) # this is true

#GAB Why? Because at some point of your net0 has a branch with 0, not 0.0. When processed by your function, it becomes 0.0 and thus different from the original string. I corrected these errors in net00 and these pass the test.
#GAB Curious about the === operator? See here https://riptutorial.com/julia-lang/example/25048/using----------and-isequal

"""
Writing the net to a file
netHyLaF = hybridlambdaformat(netScaledJ)
open("HL_3", "w") do io
	write(io, netHyLaF)
end
"""

#GAB For blocks being commented, the """ comes in handy:
"""
cd /home/carlos/hybrid-Lambda-0.6.2-beta/src
../src/hybrid-Lambda -spcu ../src/NetworkNameForHLFormat -dot -label -o NameforHyLambda > FilenameToStoreOutput.txt 2>&1
../src/hybrid-Lambda -spcu ../src/HL_3 -dot -label -o Net3 > FilenameToStoreOutput.txt 2>&1
"""
