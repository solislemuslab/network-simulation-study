#!/usr/bin/bash

Rscript Rprueba1.r
#julia JuliaPrueba.jl > JuliaHibrid2
#xdg-open NetworkjuliaR2.pdf
#xdg-open Rplots.pdf
#cp /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/JuliaHibrid2 /home/acosta/Descargas/hybrid-lambda-v0.6.1-beta/src/JuliaHibrid2
#cd /home/acosta/Descargas/hybrid-lambda-v0.6.1-beta/src
#../src/hybrid-Lambda -spcu ../src/JuliaHibrid2 -dot -label -o hybridLambdaPlot
#xdg-open hybridLambdaPlot.pdf
#cp /home/acosta/Descargas/hybrid-lambda-v0.6.1-beta/src/hybridLambdaPlot.pdf /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/
# the next script no longer needs to redirect any output, it closes silent 
julia JuliaPrueba.jl

# iterate over the content of nets_not_level1.hybridlambda and
# run hybrid-Lambda on each line in order to check issues
# with ultrametricity
while read line
do
    let LNUM+=1
    echo "Processing network in line $LNUM:"
    hybrid-Lambda -spcu $line -num 3 -seed 2 -o level2_2
done < nets_not_level1.hybridlambda
