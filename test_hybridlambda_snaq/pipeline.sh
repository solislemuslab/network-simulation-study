#!/usr/bin/bash

### Simulate gene trees from a small, level-1 network known to work with snaq
# the original one seems to be not ultrametric:
printf "Original network\n"
cat smallnet_original

# on the other hand, the modified hybrid 'major edge' do not show the warning
# but the branch length was manually modified in order to find the value which
# stops the warning, which happened to be 4.0:
printf "\nDiff between networks\n"
diff smallnet_original smallnet_ultrametric

### simulate one gene tree with hybrid-Lambda using the appropriate parameters
# for the not ultrametric version of the network
#printf "\nSimulation with the original network\n"
#hybrid-Lambdav0.6.2-pathdiff -spcu smallnet_original -num 1 -seed 1985 -o smallnet

# simulate with the manually-calibrated one
printf "\nSimulation with the manually-calibrated network\n"
hybrid-Lambdav0.6.2-pathdiff -spcu smallnet_ultrametric -num 3000 -seed 1985 -o smallnet_ultrametric

# simulate an alternative with different branch lengths and different major hybrid edge
#printf "\nSimulation with the switched major hybrid edge\n"
#hybrid-Lambdav0.6.2-pathdiff -spcu smallnet_another_major -num 1 -seed 1985 -o smallnet

# backup the original results
cp smallnet_ultrametric_coal_unit smallnet_ultrametric_coal_unit.bak

# replace the numeral with the string 'hyb'
sed -i 's/H1#0.5/H1hyb/g' smallnet_ultrametric_coal_unit

# run julia
julia estimate_net_from_gts.jl 
