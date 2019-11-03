all:

init:
	mkdir inputs
	mkdir graphs
	mkdir outputs
	mkdir errors
graph:
	python graph.py config.txt
submit:
	python submit.py config.txt submit
test: config_test.txt
	python submit.py config_test.txt

clean:
	rm outputs/*
	rm errors/*
report:
	grep "non-zero\|total time\|Degree Histogram" outputs/* | tr \,\:\(\) "\ " | tr "abcdefghijklmnopqrstuvwxyz/HD" "\ " > report.txt
