/*
 * main.h
 *
 *  Created on: May 8, 2018
 *      Author: Guang 
 */

#ifndef MAIN_H_
#define MAIN_H_

/* includes*/
#include <omp.h>
#include <iostream>
#include <climits>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <fstream>
#include <string.h>
#include <vector>
#include <assert.h>
#include <math.h>
#include <vector>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <functional>
using namespace std;
bool inter = true;
const vector<vector<bool>> setBB= {{false},{false}};
//  maxSteps,fct,ict,cct,gen,seed
const vector<vector<int>> setII= {{INT_MAX,2,100,10,0}, {INT_MAX,2,100,1,0}};
const vector<vector<double>> setDD = {{3.6, 1.0,0.5},{3.6, 1.0,0.5}};


/*problem and assignment*/
char* fileName;
int numCs;
int numVs;
int numV1;
int numC1;
int numCc;
vector<int>* clauses;
vector<int>* posC;
vector<int>* negC;
vector<int> clauseT;
uniform_int_distribution<int> distribution(0,INT_MAX);
/*methods*/

void readFile(const char* fileName);
void memAllocate(string buff);
void parseLine(string line, int index);
void initialAssignment();
void printVariables();
void printClauses();
void debugProblem();
int getFlipCandidate_wa(int cIndex);
int getFlipCandidate_max(int cIndex);
void printVector(vector<int>& vec);
void printUsage();

template<class T>
class Process{
	uniform_int_distribution<int> distribution;
	int seed;
	int* numP;
	int numUnsat = 0;
	vector<int>* unsat;
	double* probs;
	bool* assign;
	double* lookUpTable;
	int* tabuS;
	T generator;

	/*settings*/
	bool tabu_flag;
	/*option values*/


	int maxSteps;
	int fct;
	int cct;
	double cb;
	double eps;
	double lct;

public:
	Process(const vector<bool>& setB, const vector<int>& setI,const vector<double>& setD);
	void optimal();
	void parseOptions(const vector<bool>& setB, const vector<int>& setI,const vector<double>& setD);
	void printOptions();
	double LookUpTable_exp(int literal);
	double LookUpTable_poly(int literal);
	double func_exp(int literal,int partition);
	double func_poly(int literal,int partition);
	double (Process::*lookUp)(int)  = NULL;
	void initLookUpTable_exp(void);
	void initLookUpTable_poly(void);
	int getFlipLiteral(int cIndex, int partition);
	void flip(int literal);
	void flipO(int literal, int partition);
	void biasAssignment();
	void randomBiasAssignment();
	void randomAssignment();
	void setAssignment();
	void printAssignment();
	void printUnsatCs();
	void printNumP();
	void test();
	void testLine(string line);
	void debugAssign();
	int computeBreakScore(int literal,int partition);
	int (Process::*randINT)(void) = NULL;
	int randI();
	int randI2();
	void solvePart(int index);
	void push(int cIndex);
	void solve();
	void flipS(int literal);
};
#endif /* MAIN_H_ */
