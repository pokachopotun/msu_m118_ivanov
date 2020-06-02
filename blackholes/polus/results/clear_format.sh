#!/bin/bash

folder=$1

sed -i 's/.short:.*:It took me://g' $folder/time.txt
sed -i 's/.short:.*:bh_count//g' $folder/bh_count.txt
sed -i 's/.short:.*:bh_filtered//g' $folder/bh_filtered.txt
sed -i 's/outputs\///g' $folder/*.txt
