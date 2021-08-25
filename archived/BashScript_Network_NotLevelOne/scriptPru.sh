Pat1=$(pwd)
R="/RNetworks"
RJ="/RJuliaNetworks"
JH="/JuliaForHybLam"
storage="/StorageHyL"
GeneT="/GeneTrees"
UltHy="/UltraMetric/Ult_JuliaForHybLam"
UltJu="/UltraMetric/Ult_RJuliaNetworks"
UltR="/UltraMetric/Ult_RNetworks"
#
delR="$Pat1$R"
delJ="$Pat1$RJ"
delJH="$Pat1$JH"
delstorage="$Pat1$storage"
delGeneT="$Pat1$GeneT"
delUltHy="$Pat1$UltHy"
delUltJu="$Pat1$UltJu"
delUltR="$Pat1$UltR"

#
cd $delR
rm *
cd $delJ
rm *
cd $delJH
rm *
cd $delstorage
rm *
cd $delGeneT
rm *
cd $delUltHy
rm *
cd $delUltJu
rm *
cd $delUltR
rm *

cd $Pat1

while getopts t:s:n:d: flag
do
    case "${flag}" in
        t) Tree=${OPTARG};;
        s) Netsim1=${OPTARG};;
		n) num=${OPTARG};;
		d) seed=${OPTARG};;
    esac
done

echo "------------------------------------------------------------------------------"
echo "Tree: $Tree";
echo "NetSim: $Netsim1";

Rscript Rprueba1.r $Tree $Netsim1

julia JuliaPrueba.jl
#
cp /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/JuliaForHybLam/* /home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src
#
cd $delJH
a=`ls`
ini="../src/hybrid-Lambda -spcu ../src/"
subf="HybUlt"
mid=" -dot -label -o "
mid2=" > outerr"
fin=".txt 2>&1"
#
cd /home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src
for entry in $a; do
    top=${#entry}
	HybFile2=${entry:11:top}
	HybFile="$subf$HybFile2"
	res="$ini$entry$mid$subf$HybFile$mid2$HybFile$fin"
    eval $res
done
#####
del="`ls HybUlt*` `ls outerrH*` `ls JulHybrLamb*`"
patOr="mv /home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src/"
patDes=" /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/StorageHyL/"
for f in $del
do
  action="$patOr$f$patDes$f"
  eval $action
done
cd /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne
Rscript Rpru2.r
#
cp /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/UltraMetric/Ult_JuliaForHybLam/* /home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src


#
cd /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/UltraMetric/Ult_JuliaForHybLam
#
Ult=`ls`
ini="../src/hybrid-Lambda -spcu ../src/"
subf="GeneTree"
#$num
#$seed
mid=" -num $num -seed $seed -o "
#../src/hybrid-Lambda -spcu ../src/JuliaHibrid2 -num 500 -seed 2 -o example3
ro="_coal_unit"
#
patOr="mv /home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src/"
patDes=" /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/GeneTrees/"
#
cd /home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src
for entry in $Ult; do
    top=${#entry}
	HybFile2=${entry:11:top}
	res="$ini$entry$mid$subf$HybFile2"
	action="$patOr$subf$HybFile2$ro$patDes$subf$HybFile2$ro"
    eval $res
	echo "waitfora moment"
	eval $action
done
Pat1=$(pwd)
#cfSnaq

