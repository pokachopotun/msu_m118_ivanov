all: local

local: main.cpp
	g++ -std=c++11 main.cpp -o local

test: local
	./local ../testdata/chinese_article.graph -1 
