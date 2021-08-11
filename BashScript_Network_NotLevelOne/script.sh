Rscript Rprueba1.r
julia JuliaPrueba.jl > JuliaHibrid2
xdg-open NetworkjuliaR2.pdf
#xdg-open Rplots.pdf
cp /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/JuliaHibrid2 /home/acosta/Descargas/hybrid-lambda-v0.6.1-beta/src/JuliaHibrid2
cd /home/acosta/Descargas/hybrid-lambda-v0.6.1-beta/src
../src/hybrid-Lambda -spcu ../src/JuliaHibrid2 -dot -label -o hybridLambdaPlot > /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/outerr.txt 2>&1
xdg-open hybridLambdaPlot.pdf
cp /home/acosta/Descargas/hybrid-lambda-v0.6.1-beta/src/hybridLambdaPlot.pdf /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/