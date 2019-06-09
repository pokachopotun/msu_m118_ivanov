#!/bin/bash

colors=( Red Cyan Blue LightBlue Purple Yellow Lime Magenta Gray Black Orange Maroon Green Olive )

filetype=$1;
email=$2;

files=$(find . -name "*${filetype}" -type f)
echo ${filetype}

function change {
	
	IFS='/' read -r -a old <<< "$1"
	IFS='/' read -r -a new <<< "$2"

	oldLen=${#old[@]}
	newLen=${#new[@]}
	
	comLen=0;
	for (( i=0; i < $newLen && i < $oldLen; i++ ));
	do
		comLen=$i
		if [[ ${old[$i]} != ${new[$i]} ]]
		then
			break
		fi
	done

	for (( i=$oldLen - 1; i > $comLen; i-- ));
	do
		echo "</ul>"
		echo "</li>"
	done	

	for (( i=$comLen; i < $newLen - 1; i++ ));
	do
		style=$( getStyle $i )
		echo "<li $style> ${new[$i]}"
		echo "<ul>"
	done	
	
	echo "<li $( getStyle $(( $newLen - 1 )) )> ${new[$newLen-1]} </li>"
}

function getStyle {
	depth=$1
	cnt=${#colors[@]}
	cid=$(( $depth % $cnt ))
	echo "style=\"color:${colors[$cid]};\""
}

function close {
	IFS='/' read -r -a path <<< "$1"

	len=${#path[@]}
	for (( i=0; i < len; i++));
	do
		echo "</ul>"
		echo "</li>"
	done
}



echo "<ul>"
curPath="."
for path in $files
do
	change $curPath $path
	curPath=$path
	#echo $curPath
	#echo ${len} ${path[@]}
done

close curPath

echo "</ul>"

