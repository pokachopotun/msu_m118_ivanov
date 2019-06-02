#!/bin/bash

case "$1" in
[1] ) mkdir NEW; echo "1) NEW dir created";;
[2] ) echo "2) Current directory is: "; pwd;;
[3] ) echo "3) Version info: "; cat /proc/version;;
[4] ) echo "4) Print all txts' content: "; cat *.txt;;
[5] ) echo "5) Print directories: "; ls */ -d;;
[6] ) echo "6) Mounted devices and emory usage: "; df -h;;
[7] ) echo "7) Environment variables: "; printenv;;
* ) echo "There is no such command";;
esac
