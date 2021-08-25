while getopts n:m:s:t: flag
do
    case "${flag}" in
        n) num=${OPTARG};;
        m) mean=${OPTARG};;
        s) sd=${OPTARG};;
		t) text=${OPTARG};;
    esac
done

echo "------------------------------------------------------------------------------"
echo "Numero de observaciones: $num";
echo "mean: $mean";
echo "sd: $sd";
echo "texto: $text";

Rscript Try1.r $num $mean $sd $text