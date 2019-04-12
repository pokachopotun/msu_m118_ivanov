import os
import sys

if __name__ == "__main__":
	filePath = sys.argv[1];
	cnt = int( sys.argv[2] );
	print(cnt)
	with open(filePath, 'rb' ) as file:
		content = file.read();
	line = ""
	c = 0
	for i in range(len(content)):
		if i % cnt == 0:
			print(line)
			line = ""
			c+=1
		if int(content[i]) == 0:
			line += "0"
		else:
			line += "1"
	print(line)
	print(c)
