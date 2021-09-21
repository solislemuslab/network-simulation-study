Path0=$(pwd)
cd ..
Pat1=$(pwd)

mkdir -p LastRunnings
n=`ls LastRunnings | wc -l`
n1=$(($n+1))
a=`ls`
cu="Output"
for entry in $a; do
	if [ $entry == $cu ]
	then
	 mkdir -p "LastRunnings/Output$n1"
	 mv Output/* "LastRunnings/Output$n1"
	 rmdir Output
	fi
done

#Path to src HybridLambda
PatSrcHyLa="/home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src"

#Paths to save temporarily the results
R="Output/RNetworks"
RJ="Output/RJuliaNetworks"
JHL="Output/JuliaForHybLam"
TstoHL="Output/TemporaryHyLambdaStorage"
CF="Output/SnaqOut/CFactors"
GT="Output/SnaqOut/GeneTrees"
NET="Output/SnaqOut/SnaqNet"
patGT="$Pat1/$GT"
#RNetWar="Output/SnaqOut/RNetWarnings"


mkdir -p {$R,$RJ,$JHL,$TstoHL,$CF,$GT,$NET}
pathR=$Pat1/$R
pathRJ=$Pat1/$RJ
pathJHL=$Pat1/$JHL
patStoHL=$Pat1/$TstoHL
#patRNetWar=$Pat1/$RNetWar 


#paths of ultrametric Networks
UlJHL="Output/UltraMetric/Ult_JuliaForHybLam"
UltRJ="Output/UltraMetric/Ult_RJuliaNetworks"
UltR="Output/UltraMetric/Ult_RNetworks"
UltGT="Output/UltraMetric/GeneTrees"
UltHLStor="Output/UltraMetric/UltHyLambdaStorage"
mkdir -p {$UlJHL,$UltRJ,$UltR,$UltGT,$UltHLStor}


PatUltR=$Pat1/$UltR
PatUltRJ=$Pat1/$UltRJ
PatUlJHL=$Pat1/$UlJHL
PatUltGT=$Pat1/$UltGT
PatUltHLStor=$Pat1/$UltHLStor


#paths of Non ultrametric Networks
NuJHL="Output/NonUltraMetric/NonUlt_JuliaForHybLam"
NuRJ="Output/NonUltraMetric/NonUlt_RJuliaNetworks"
NuR="Output/NonUltraMetric/NonUlt_RNetworks"
NuHLSt="Output/NonUltraMetric/NonUltHyLambdaStorage"
mkdir -p {$NuJHL,$NuRJ,$NuR,$NuHLSt}

PatNuJHL=$Pat1/$NuJHL
PatNuRJ=$Pat1/$NuRJ
PatNuR=$Pat1/$NuR
PatNuHLSt=$Pat1/$NuHLSt





cd $Path0

while getopts t:s:n:d: flag
do
    case "${flag}" in
        t) Tree=${OPTARG};;
        s) Netsim1=${OPTARG};;
		n) num=${OPTARG};;
		d) seed=${OPTARG};;
    esac
done

RANDOM=$seed
HyLaSeed="$((1 + $RANDOM % 1000))"
Rseed="$((1 + $RANDOM % 1000))"


Rscript SimSpeciesNetworks.r "$Tree" "$Netsim1" "$pathR" "$Rseed"
julia Net_HybridLamb_Format.jl "$pathR" "$pathRJ" "$pathJHL" 


#Copy SpeciesNetworks to hybrid lambda SRC
ini="cp "
med="/* "
la="$ini$pathJHL$med$PatSrcHyLa"
eval $la



cd $pathJHL
a=`ls`
ini="../src/hybrid-Lambda -spcu ../src/"
subf="HybUlt"
mid=" -dot -label -o "
mid2=" > outerr"
fin=".txt 2>&1"
#
cd $PatSrcHyLa
for entry in $a; do
    top=${#entry}
	HybFile2=${entry:11:top}
	HybFile="$subf$HybFile2"
	res="$ini$entry$mid$subf$HybFile$mid2$HybFile$fin"
    eval $res
done

#List of all the files that Hybrid lambda may create
del="`ls HybUlt*` `ls outerrH*` `ls JulHybrLamb*`"
patOr="mv ""$PatSrcHyLa""/"
patDes=" ""$patStoHL""/"
for f in $del
do
  action="$patOr$f$patDes$f"
  eval $action
done





#cd /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne
cd $Path0


Rscript UltrametricEvaluation.r $patStoHL $pathR $pathRJ $pathJHL $PatUltR $PatUltRJ $PatUlJHL $PatUltHLStor $PatNuR $PatNuRJ $PatNuJHL $PatNuHLSt




#How to delet the temporarily  files
rm -r $pathR
rm -r $pathRJ
rm -r $pathJHL
rm -r $patStoHL




#Copy Ultrametric SpeciesNetworks to hybrid lambda SRC
ini="cp "
med="/* "
la="$ini$PatUlJHL$med$PatSrcHyLa"
eval $la


#
cd $PatUlJHL
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

ini1="mv "
Fins="/"
#patOr="mv /home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src/"
patOr="$ini1$PatSrcHyLa$Fins"
ini2=" "
patDes="$ini2$PatUltGT$Fins"
#
cd $PatSrcHyLa
for entry in $Ult; do
    top=${#entry}
	HybFile2=${entry:11:top}
	res="$ini$entry$mid$subf$HybFile2"
	action="$patOr$subf$HybFile2$ro$patDes$subf$HybFile2$ro"
    eval $res
	echo "waitfora moment"
	eval $action
done

#Copy Ultrametric SpeciesNetworks to hybrid lambda SRC
ini="cp "
med="/* "
la="$ini$PatUltGT$med$patGT"
eval $la
#

