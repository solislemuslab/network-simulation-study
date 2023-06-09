
path_ini=`pwd`

cd ../data

tar_name="files.tar.gz"
index_name="chtc"

# Estoy dentro de data
dir_store="../Files_chtc/$index_name"
ini=-1
for f in *; do
    ini=`expr $ini + 1`
    dir_copy="$dir_store$ini"
    mkdir -p "$dir_copy"
    #file= eval "ls $f"
    cd $f
    fun="tar -czvf $tar_name `ls`"
    eval $fun
    cd ..
    eval $"cp $f/$tar_name $dir_copy"
    eval $"rm $f/$tar_name"
done

eval $"cd $path_ini"










