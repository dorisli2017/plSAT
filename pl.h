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
bool sat0 = false;
bool sat2 = false;
bool sat = false;
bool* assignG;
//  maxSteps,fct,ict,cct,gen,seed
const vector<vector<int>> setII= {{INT_MAX,2,100,10,0}, {INT_MAX,2,100,10,0}};
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
int maxL = 0;
/*methods*/

void readFile(const char* fileName);
void memAllocate(string buff);
void parseLine(string line, int index);
void initialAssignment();
void printVariables();
void printClauses();
void debugProblem();
void printVector(vector<int>& vec);
void printUsage();
void testLine(string line);
void testPart(int partition, bool* assignG);

template<class T>
class Process{
	uniform_int_distribution<int> distribution;
	int seed;
	int* numP;
	vector<int> unsat;
	double* probs;
	bool* assign;
	double* lookUpTable;
	int* tabuS;
	T generator;
	double noise = 2.0;
	int *breaks;
	int *critVar;

	/*settings*/
	/*option values*/
	int fct = 1;
	int cct;
	double cb;
	double eps;
	double lct;
	long long int flipCount = 0;
	vector<int> clauseQ;

public:
	Process();
	void optimal();
	void printOptions();
	double LookUpTable_exp(int literal);
	double LookUpTable_poly(int literal);
	double (Process::*lookUp)(int)  = NULL;
	int (Process::*computeBreak)(int) = NULL;
	void initLookUpTable_exp(void);
	void initLookUpTable_poly(void);
	int (Process::*getFlipLiteral)(int cIndex, int partition);
	int getFlipLiteral3(int cIndex, int partition);
	int getFlipLiteral57(int cIndex, int partition);
	void (Process::*flip)(int literal, int partition);
	void flip3(int literal, int partition);
	void flip57(int literal, int partition);
	void (Process::*initAssignment)(int) = NULL;
	void biasAssignment(int partition);
	void randomAssignment(int partition);
	void (Process::*setAssignment)(int partition) = NULL;
	void setAssignment3(int partition);
	void setAssignment57(int partition);
	void printAssignment();
	void printUnsatCs();
	void printNumP();
	void test();
	void debugAssign();
	int computeBreakScore(int literal);
	int computeBreakScore0(int literal);
	int computeBreakScore2(int literal);
	int (Process::*randINT)(void) = NULL;
	int randI();
	int randI2();
	void solvePart(int index);
	void push(int cIndex);
	void solve();
};
#endif /* MAIN_H_ */
