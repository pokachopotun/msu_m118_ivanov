#!/bin/bash

sed -i 's/.short:.:It took me://g' time.txt
sed -i 's/.short:.:bh_count//g' bh_count.txt
sed -i 's/.short:.:bh_filtered//g' bh_filtered.txt
sed -i 's/outputs\///g' *.txt
