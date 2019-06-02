#!/bin/bash

filename=$1

echo $filename

while IFS= read -r line
do
	IFS=';' read -r -a params <<< "$line"
	newfile=${params[0]}
	user=${params[1]}
	group=${params[2]}
	access=${params[3]}
	touch ${newfile}
	chown ${user}:${group} ${newfile}
	if [ $? != 0 ]
		then 
			rm ${newfile}
	fi
#	echo "$newfile"
#	echo "$user"
#	echo "$group"
#	echo "$access"
done < "$filename"
