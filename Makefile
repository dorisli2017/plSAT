all:
	g++ -O3 -fopenmp pl.cpp -std=c++11 -o pl
clean:
	rm pl
