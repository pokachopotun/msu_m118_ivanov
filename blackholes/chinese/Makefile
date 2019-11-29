all: local

local: main.cpp
	g++ -O3 -std=c++11 main.cpp -o local

test: local
	./local ../testdata/chinese_article.graph -1 
