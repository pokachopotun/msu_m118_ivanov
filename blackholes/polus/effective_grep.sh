#!/bin/bash

folder=$1
results=$2
rm -r $folder/*.short

total=$(ls $folder/ | wc -l)
echo $total
cnt=0
for file in $folder/*
do
    cnt=$(( $cnt + 1 ))
    echo "$cnt / $total"
    newfile=$file.short
    touch $newfile
    tac $file | grep -m 1 "bh_count" >> $newfile
    tac $file | grep -m 1 "bh_filtered" >> $newfile
    tac $file | grep -m 1 "It took me" >> $newfile
done
