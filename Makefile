all:
	g++ -O3 -fopenmp pl.cpp -std=c++11 -o pl2
clean:
	rm pl
# pl1 two threads no tabu (max 1)
#pl2 two thread tabu (max 1)
