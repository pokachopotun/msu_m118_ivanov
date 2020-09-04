first=$1
second=$2
#third=$3
output=$3
mkdir $output
for file in $first/*.txt
do
    name=$(basename $file)
    cat $file $second/$name > $output/$name
done
