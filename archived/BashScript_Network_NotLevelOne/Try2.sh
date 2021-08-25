while getopts u:a:f: flag
do
    case "${flag}" in
        u) username=${OPTARG};;
        a) age=${OPTARG};;
        f) fullname=${OPTARG};;
    esac
done

echo "Username: $username";
echo "Age: $age";
echo "Full Name: $fullname";

Rscript Try1.r $username $age $fullname