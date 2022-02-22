#Add General description as the readme file


#Get the Path where there are the scrips | PathScripts= PATH_SCRIPTS = PATH_SCRIPTS  PARENT_PTH_SCRIPTS 
PATH_SCRIPTS=$(pwd)

# Get parent directory of the "path of the scrips" (In this path we will create all the directories nedded to save the results) | ParentPathScripst = PATH_SCRIPTS
cd ..
PARENT_PTH_SCRIPTS=$(pwd)


#Path to src HybridLambda (Fixed Path)
PTH_HYLA="/home/carlos/hybrid-Lambda-0.6.2-beta/src"

#Pat for MaxCut (Fixed Path)
PTH_QMCT="/home/carlos/phyloNet/PhyloNetworks.jl.wiki/data_results/scripts"



#Path to save the rNetworks with warnings TEMP_PTH_R_WARN=PTH_WARN_R_OUT
R_WARN="Output/RNetworksWarnings"
PTH_WARN_R_OUT="$PARENT_PTH_SCRIPTS/$R_WARN"
mkdir -p "$R_WARN"


#Path to save The Snaq input (Our final results in this part)
NET="Output/snaq_input"
mkdir -p "$NET"
PTH_SNAQ_INPUT=$PARENT_PTH_SCRIPTS/$NET




#Paths to save temporarily the results
R="Output/RNetworks"
RJ="Output/RJuliaNetworks"
JHL="Output/JuliaForHybLam"
TstoHL="Output/TemporaryHyLambdaStorage"
mkdir -p {$R,$RJ,$JHL,$TstoHL,$NET}
TEMP_PTH_R=$PARENT_PTH_SCRIPTS/$R
TEMP_PTH_RJ=$PARENT_PTH_SCRIPTS/$RJ
TEMP_PTH_JHL=$PARENT_PTH_SCRIPTS/$JHL
TEMP_PTH_STORE_HL=$PARENT_PTH_SCRIPTS/$TstoHL


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
PATH_ULT_R=$PARENT_PTH_SCRIPTS/$UltR
PATH_ULT_RJ=$PARENT_PTH_SCRIPTS/$UltRJ
PATH_ULT_JHL=$PARENT_PTH_SCRIPTS/$UlJHL
PATH_ULT_GT=$PARENT_PTH_SCRIPTS/$UltGT
PATH_ULT_HYLA_STORAGE=$PARENT_PTH_SCRIPTS/$UltHLStor
PATH_ULT_README=$PARENT_PTH_SCRIPTS/$UltReadme
PATH_ULT_C_FACTOR=$PARENT_PTH_SCRIPTS/$UltCFact
PATH_ULT_QMCT=$PARENT_PTH_SCRIPTS/$UltQMCT



#paths of Non ultrametric Networks
NuJHL="Output/NonUltraMetric/NonUlt_JuliaForHybLam"
NuRJ="Output/NonUltraMetric/NonUlt_RJuliaNetworks"
NuR="Output/NonUltraMetric/NonUlt_RNetworks"
NuHLSt="Output/NonUltraMetric/NonUltHyLambdaStorage"
mkdir -p {$NuJHL,$NuRJ,$NuR,$NuHLSt}
PATH_NOULT_JHL=$PARENT_PTH_SCRIPTS/$NuJHL
PATH_NOULT_RJ=$PARENT_PTH_SCRIPTS/$NuRJ
PATH_NOULT_R=$PARENT_PTH_SCRIPTS/$NuR
PATH_NOULT_HYLA_STOR=$PARENT_PTH_SCRIPTS/$NuHLSt




########
#start the process

#Arguments used for run the Species_Gene_Networks.sh
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

RSEED="$((1 + $RANDOM % 1000))"



#Go to the inicial path (scripth path)
cd $PATH_SCRIPTS


#TimeSECONDS=0
Rscript SimSpeciesNetworks.r "$Tree" "$Netsim1" "$TEMP_PTH_R" "$RSEED" "$PTH_WARN_R_OUT"
#TimeRduration=$SECONDS




#convert the networks to hybrid lambda format
#TimeSECONDS=0
julia Net_HybridLamb_Format.jl "$TEMP_PTH_R" "$TEMP_PTH_RJ" "$TEMP_PTH_JHL" 
#TimeJuduration=$SECONDS




#TimeSECONDS=0
#Copy SpeciesNetworks to hybrid lambda SRC
#cp /home/carlos/GitProjects/PhyloNetworks/network-simulation-study/snaq_level2/Output/JuliaForHybLam/* /home/carlos/hybrid-Lambda-0.6.2-beta/src
ini="cp "
med="/* "
la="$ini$TEMP_PTH_JHL$med$PTH_HYLA"
eval $la


#Run hybrid lambda for each network to know if is UTRAMETRIC or not
#../src/hybrid-Lambda -spcu ../src/NetworkNameForHLFormat -dot -label -o NameforHyLambda > FilenameToStoreOutput.txt 2>&1
cd $TEMP_PTH_JHL
NAM_NETWORK_HL_FORMAT=`ls`

HYLA_CODE_INIT="../src/hybrid-Lambda -spcu ../src/"
SUBNAME_HYLA="HybUlt"
HYLA_CODE_NOMINATE=" -dot -label -o "
NOMINATE_ULT_OUT=" > outerr"
HYLA_CODE_FIN=".txt 2>&1"

cd $PTH_HYLA
for NAM_NET in $NAM_NETWORK_HL_FORMAT; do
    top=${#NAM_NET}
	HybFile2=${NAM_NET:11:top}
	NETWORK_ID="$SUBNAME_HYLA$HybFile2"
	res="$HYLA_CODE_INIT$NAM_NET$HYLA_CODE_NOMINATE$SUBNAME_HYLA$NETWORK_ID$NOMINATE_ULT_OUT$NETWORK_ID$HYLA_CODE_FIN"
    eval $res
done


#List of all the files that Hybrid lambda may create to remove them so that hybrid-lambda-src is as at the beginning
del="`ls HybUlt*` `ls outerrH*` `ls JulHybrLamb*`"

patOr="mv ""$PTH_HYLA""/"
patDes=" ""$TEMP_PTH_STORE_HL""/"
for f in $del
do
  action="$patOr$f$patDes$f"
  eval $action
done

#TimeHLduration=$SECONDS



#cd /home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne
cd $PATH_SCRIPTS
#TimeSECONDS=0
Rscript UltrametricEvaluation.r $TEMP_PTH_STORE_HL $TEMP_PTH_R $TEMP_PTH_RJ $TEMP_PTH_JHL $PATH_ULT_R $PATH_ULT_RJ $PATH_ULT_JHL $PATH_ULT_HYLA_STORAGE $PATH_NOULT_R $PATH_NOULT_RJ $PATH_NOULT_JHL $PATH_NOULT_HYLA_STOR
#TimeRUltDuration=$SECONDS





#Delet the temporarily  files
rm -r $TEMP_PTH_R
rm -r $TEMP_PTH_RJ
rm -r $TEMP_PTH_JHL
rm -r $TEMP_PTH_STORE_HL





#Get Gene Trees
#TimeSECONDS=0
#Copy Ultrametric SpeciesNetworks to hybrid lambda SRC to get GeneTrees
ini="cp "
med="/* "
la="$ini$PATH_ULT_JHL$med$PTH_HYLA"
eval $la


cd $PATH_ULT_JHL
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
patOr="$ini1$PTH_HYLA$Fins"
ini2=" "
patDes="$ini2$PATH_ULT_GT$Fins"

#Modificar luego
patDes2="$ini2$PATH_ULT_README$Fins"
lo="ju.txt"
Rea="readme"
filis=".txt"
##
##
cd $PTH_WARN_R_OUT
nWar=`ls | wc -l`
cd $PATH_ULT_R
nUlt=`ls | wc -l`
cd $PATH_NOULT_R
nNotUlt=`ls | wc -l`

#Modificar luego

cd $PTH_HYLA
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
	echo $Tree" "$Netsim1"  NumUlt = "$nUlt"  NumNotUlt = "$nNotUlt";  GlobalSeed: "$seed";  RSEED: "$RSEED";  GTSeed: "$HyLaSeed>ju.txt
	action2="$patOr$lo$patDes2$Rea$HybFile2$filis"
	eval $action2
done

#TimeGTDuration=$SECONDS















stringZ="$Tree"
n1=`expr index "$stringZ" =`
lo=";"
n2=`expr index "$stringZ" $lo`

RnumSim=${stringZ:$n1:$n2-$n1-1}
#####
#Timeduration=$SECONDS


echo "-----------------------------------------------------"
echo "NumSIm: $RnumSim"
echo "NumR_WARNing: $nWar"
echo "NumRUlt: $nUlt"
echo "NumRNotUlt: $nNotUlt"
echo "-----------------------------------------------------"

#TIME of the process
#Timeecho "R$(($Rduration / 60)) minutes and $(($Rduration % 60)) seconds elapsed."
#Timeecho "Ju$(($Juduration / 60)) minutes and $(($Juduration % 60)) seconds elapsed."
#Timeecho "HL$(($HLduration / 60)) minutes and $(($HLduration % 60)) seconds elapsed."
#Timeecho "RUlt$(($RUltDuration / 60)) minutes and $(($RUltDuration % 60)) seconds elapsed."
#Timeecho "GT$(($GTDuration / 60)) minutes and $(($GTDuration % 60)) seconds elapsed."


#Get all the the files nescessary for run Snaq (Concordance Factors, Consensous Tree (obtained with QMCT)), and ohres complementary files.
#Then storage it in a suitable path
cd $PATH_SCRIPTS
julia SnaqHybridzationPreparing.jl "$PATH_ULT_C_FACTOR" "$PATH_ULT_GT" "$PTH_QMCT" "$PATH_ULT_QMCT" "$PATH_ULT_R" "$PATH_ULT_JHL" "$PTH_SNAQ_INPUT" "$PATH_ULT_README"



#Nominate the final output with the parameters chosen
a=`date +'%m-%d-%Y-%H:%M:%S'`
MainFileNam1=$Tree"_"$Netsim1"_GTnum="$num"_Date:"$a
T1=`echo $MainFileNam1 | tr ";""="":""," _-`
one="${T1/c(/}"
MainFileNam="${one/)/}"


cd $PARENT_PTH_SCRIPTS
mv Output $MainFileNam

