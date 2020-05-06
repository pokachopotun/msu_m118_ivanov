#!/bin/bash

folder=$1
config=$1.conf

mkdir $folder
cd $folder
mkdir outputs
mkdir errors
cd ..

runsh=$folder/run.sh

python run.py $config > $runsh
chmod +x $runsh
