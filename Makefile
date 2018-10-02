all:
	g++ -fopenmp -O3 a1.cpp -std=c++11 -o a1
clean:
	rm pl
# pl1 two threads no tabu (max 1)
#pl2 two thread tabu (max 1)
