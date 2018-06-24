all:
	g++ -O3 -fopenmp l.cpp -std=c++11 -o l
clean:
	rm l
