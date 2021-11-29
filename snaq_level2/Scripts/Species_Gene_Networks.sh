
# do some work


Path0=$(pwd)
cd ..
Pat1=$(pwd)



#Path to src HybridLambda
#PatSrcHyLa="/home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src"
PatSrcHyLa="/home/carlos/hybrid-Lambda-0.6.2-beta/src"

#Pat for MaxCut
PatMaxCu="/home/carlos/phyloNet/PhyloNetworks.jl.wiki/data_results/scripts"
#FixedPath


#Path to save the rNetworks with warnings
Rwarn="Output/RNetworksWarnings"
patRwarn="$Pat1/$Rwarn"


#Paths to save temporarily the results
R="Output/RNetworks"
RJ="Output/RJuliaNetworks"
JHL="Output/JuliaForHybLam"
TstoHL="Output/TemporaryHyLambdaStorage"
#CF="Output/SnaqOut/CFactors"
#GT="Output/SnaqOut/GeneTrees"
NET="Output/SnaqOut/SnaqNet"
#patGT="$Pat1/$GT"


mkdir -p {$R,$RJ,$JHL,$TstoHL,$NET,$Rwarn}
pathR=$Pat1/$R
pathRJ=$Pat1/$RJ
pathJHL=$Pat1/$JHL
patStoHL=$Pat1/$TstoHL
patNET=$Pat1/$NET
#patRNetWar=$Pat1/$RNetWar 


#paths of ultrametric Networks
UlJHL="Output/UltraMetric/Ult_JuliaForHybLam"
UltRJ="Output/UltraMetric/Ult_RJuliaNetworks"
UltR="Output/UltraMetric/Ult_RNetworks"
UltGT="Output/UltraMetric/GeneTrees"

UltCFact="Output/UltraMetric/CFactors"
UltQMCT="Output/UltraMetric/QMCtree"

UltReadme="Output/UltraMetric/Readme"

UltHLStor="Output/UltraMetric/UltHyLambdaStorage"
mkdir -p {$UlJHL,$UltRJ,$UltR,$UltGT,$UltHLStor,$UltCFact,$UltQMCT,$UltReadme}


PatUltR=$Pat1/$UltR
PatUltRJ=$Pat1/$UltRJ
PatUlJHL=$Pat1/$UlJHL
PatUltGT=$Pat1/$UltGT
PatUltHLStor=$Pat1/$UltHLStor

PatUltReadme=$Pat1/$UltReadme

PatUltCFact=$Pat1/$UltCFact
PatUltQMCT=$Pat1/$UltQMCT

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

Rseed="$((1 + $RANDOM % 1000))"





SECONDS=0
Rscript SimSpeciesNetworks.r "$Tree" "$Netsim1" "$pathR" "$Rseed" "$patRwarn"
Rduration=$SECONDS





SECONDS=0
julia Net_HybridLamb_Format.jl "$pathR" "$pathRJ" "$pathJHL" 
Juduration=$SECONDS




SECONDS=0
#Copy SpeciesNetworks to hybrid lambda SRC
ini="cp "
med="/* "
la="$ini$pathJHL$med$PatSrcHyLa"
eval $la

#../src/hybrid-Lambda -spcu ../src/NetworkNameForHLFormat -dot -label -o NameforHyLambda > FilenameToStoreOutput.txt 2>&1

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

HLduration=$SECONDS



#cd /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne
cd $Path0
SECONDS=0
Rscript UltrametricEvaluation.r $patStoHL $pathR $pathRJ $pathJHL $PatUltR $PatUltRJ $PatUlJHL $PatUltHLStor $PatNuR $PatNuRJ $PatNuJHL $PatNuHLSt
RUltDuration=$SECONDS





#Delet the temporarily  files
rm -r $pathR
rm -r $pathRJ
rm -r $pathJHL
rm -r $patStoHL





#Get Gene Trees
SECONDS=0
#Copy Ultrametric SpeciesNetworks to hybrid lambda SRC to get GeneTrees
ini="cp "
med="/* "
la="$ini$PatUlJHL$med$PatSrcHyLa"
eval $la


cd $PatUlJHL
#
Ult=`ls`
ini="../src/hybrid-Lambda -spcu ../src/"
subf="GeneTree"
mid=" -num $num -seed "
#../src/hybrid-Lambda -spcu ../src/JuliaHibrid2 -num 500 -seed 2 -o example3
ro="_coal_unit"


ini1="mv "
Fins="/"
#patOr="mv /home/acosta/Descargas/hybrid-Lambda-0.6.2-beta/src/"
patOr="$ini1$PatSrcHyLa$Fins"
ini2=" "
patDes="$ini2$PatUltGT$Fins"

#Modificar luego
patDes2="$ini2$PatUltReadme$Fins"
lo="ju.txt"
Rea="readme"
filis=".txt"
##
##
cd $patRwarn
nWar=`ls | wc -l`
cd $PatUltR
nUlt=`ls | wc -l`
cd $PatNuR
nNotUlt=`ls | wc -l`

#Modificar luego

cd $PatSrcHyLa
for entry in $Ult; do
    HyLaSeed="$((1 + $RANDOM % 1000))"
	g=" -o "
	mid2="$mid$HyLaSeed$g"
    top=${#entry}
	HybFile2=${entry:11:top}
	res="$ini$entry$mid2$subf$HybFile2"
	action="$patOr$subf$HybFile2$ro$patDes$subf$HybFile2$ro"
    eval $res
	echo "waitfora moment"
	eval $action
	echo $Tree" "$Netsim1"  NumUlt = "$nUlt"  NumNotUlt = "$nNotUlt";  GlobalSeed: "$seed";  Rseed: "$Rseed";  GTSeed: "$HyLaSeed>ju.txt
	action2="$patOr$lo$patDes2$Rea$HybFile2$filis"
	eval $action2
done

GTDuration=$SECONDS





action2="$patOr$lo$patDes$Rea$HybFile2$filis"










stringZ="$Tree"
n1=`expr index "$stringZ" =`
lo=";"
n2=`expr index "$stringZ" $lo`

RnumSim=${stringZ:$n1:$n2-$n1-1}
#####
duration=$SECONDS


echo "-----------------------------------------------------"
echo "NumSIm: $RnumSim"
echo "NumRWarning: $nWar"
echo "NumRUlt: $nUlt"
echo "NumRNotUlt: $nNotUlt"
echo "-----------------------------------------------------"


echo "R$(($Rduration / 60)) minutes and $(($Rduration % 60)) seconds elapsed."
echo "Ju$(($Juduration / 60)) minutes and $(($Juduration % 60)) seconds elapsed."
echo "HL$(($HLduration / 60)) minutes and $(($HLduration % 60)) seconds elapsed."
echo "RUlt$(($RUltDuration / 60)) minutes and $(($RUltDuration % 60)) seconds elapsed."
echo "GT$(($GTDuration / 60)) minutes and $(($GTDuration % 60)) seconds elapsed."

cd $Path0
julia SnaqHybridzationPreparing.jl "$PatUltCFact" "$PatUltGT" "$PatMaxCu" "$PatUltQMCT" "$PatUltR" "$PatUlJHL" "$patNET" "$PatUltReadme"



#get name with the parameters chosen
a=`date +'%m-%d-%Y-%H:%M:%S'`
MainFileNam1=$Tree"_"$Netsim1"_GTnum="$num"_Date:"$a
T1=`echo $MainFileNam1 | tr ";""="":""," _-`
one="${T1/c(/}"
MainFileNam="${one/)/}"


cd $Pat1
mv Output $MainFileNam

