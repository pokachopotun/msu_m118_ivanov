#!/bin/bash

../effective_grep.sh outputs/ .

grep -n "It took me" outputs/*.short > time.txt
grep -n "bh_count" outputs/*.short > bh_count.txt
grep -n "bh_filtered" outputs/*.short > bh_filtered.txt

resdir=../results/$(basename $(pwd))
mkdir $resdir
cp *.txt $resdir
