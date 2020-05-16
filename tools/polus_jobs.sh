#!/bin/bash

tmpfile=$(mktemp)
ssh -i ~/.ssh/edu-cmc-sqi16-14 edu-cmc-sqi16-14@polus.hpc.cmc.msu.ru '( bjobs | wc -l )' &> $tmpfile && tail -n 1 $tmpfile
rm $tmpfile
