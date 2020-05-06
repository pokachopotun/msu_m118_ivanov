first=$1
second=$2
output=$3
mkdir $output
for file in $1/*.txt
do
    name=$(basename $file)
    cat $file $second/$name > $output/$name
done
