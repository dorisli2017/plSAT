/*
 * main.cpp
 *
 *  Created on: May 8, 2018
 *      Author: ping
 */

#include "pl.h"

int main(int argc, char *argv[]){
	fileName = argv[1];
	readFile(fileName);
int tid;
#pragma omp parallel num_threads(2) private(tid)
 {
	//tid = omp_get_thread_num();
	tid = 0;
	const vector<bool>setB = setBB[tid];
	const vector<int> setI =setII[tid];
	const vector<double>& setD = setDD[tid];
	switch(setI[3]){
	case 0:{
		Process<minstd_rand0> process (setB, setI,setD);
		process.optimal();
		break;
	}
	case 1:{
		Process<minstd_rand> process (setB, setI,setD);
		process.optimal();
		break;
	}
	case 2:{
		Process<mt19937> process (setB, setI,setD);
		process.optimal();
		break;
	}
	case 3:{
		Process<mt19937_64> process (setB, setI,setD);
		process.optimal();
		break;
	}
	case 4:{
		Process<ranlux24_base> process (setB, setI,setD);
		process.optimal();
		break;
	}
	case 5:{
		Process<ranlux48_base> process (setB, setI,setD);
		process.optimal();
		break;
	}
	case 6:{
		Process<ranlux24> process (setB, setI,setD);
		process.optimal();
		break;
	}
	case 7:{
		Process<ranlux48> process (setB, setI,setD);
		process.optimal();
		break;
	}
	case 8:{
		Process<knuth_b> process (setB, setI,setD);
		process.optimal();
		break;
	}
	default:{
		Process<default_random_engine> process (setB, setI,setD);
		process.optimal();
		break;
	}
	}
	//randTest();
	//debugProblem();
	//process.printOptions()
}	//process.debugAssign();
}
void debugProblem(){
	printVariables();
	printClauses();
}
template<class T>
void Process<T>::debugAssign(){
	/* Testing code**********************************/
		printOptions();
	   	printAssignment();
	   	printUnsatCs();
	   	printNumP();


	/*Testing code**********************************/

}
template<class T>
Process<T>::Process(const vector<bool>& setB, const vector<int>& setI,const vector<double>& setD){
	parseOptions(setB, setI,setD);
	distribution = uniform_int_distribution<int>(0,INT_MAX);
	//set the parameters
	   // set tabuS
	if(setI[7] == 10){
		srand(seed);
		randINT = &Process::randI2;
	}
	else {
		randINT = &Process::randI;
		generator.seed(seed);
	}
	if(tabu_flag){
		tabuS = (int*) malloc(sizeof(int) * numVs);
		for(int i = 0; i < numVs; i++){
			tabuS[i] = 0;
		}
	}
	numP = (int*) malloc(sizeof(int) * numCs);
	probs = (double*)malloc(sizeof(double) * numVs);
	assign = (bool*)malloc(sizeof(bool) * numVs);
	unsat= new vector<int>[3];
	biasAssignment();
	//set lookuptable
	switch (fct){
	case 0:initLookUpTable_poly();
			lookUp =&Process::LookUpTable_poly;
			break;
	default:initLookUpTable_exp();
			lookUp =&Process::LookUpTable_exp;
			break;
	}
}

/*parse the argument (including options and filename)
 *using getopt_long to allow GNU-style long options as well as single-character options
 */
template<class T>
void Process<T>::parseOptions(const vector<bool>& setB, const vector<int>& setI,const vector<double>& setD){
	tabu_flag = setB[0];

	maxSteps = setI[0];
	fct= setI[1];
	cct= setI[2];
	seed = setI[4];
	cb=setD[0];
	eps= setD[1];
	lct = setD[2];
}
// construct the Problem with fill information of the input file
void readFile(const char* fileName){
	//cout<< "in readFile"<<endl;
	ifstream fp;
	//fp.open("Debug/instance.cnf",std::ios::in);
	fp.open(fileName,std::ios::in);
	if(!fp.is_open()){
		perror("read file fails");
		exit(EXIT_FAILURE);
	}
	string buff;
	char head;
   	getline(fp,buff);
   	// Get the p line
   	while(!fp.eof()){
		//cout<<buff<<endl;
		//todo:parseLine
   		if(buff.empty()) continue;
		head =buff.at(0);
		if(head == 'p'){
			memAllocate(buff);
			break;
		}
	  getline(fp,buff);
	}
	if(inter){
		getline(fp,buff);
		head =buff.at(0);
		if(head == 'c'){
		}
	}// for partition file
   	// Get the clauses
   	int index = -1;
   	int line = 0;
   	while(!fp.eof() && line < numCs){
   		index++;
		getline(fp,buff);
		if(buff.empty()) continue;
		parseLine(buff, index);
		line++;
   	}
   	initialAssignment();
};
void memAllocate(string buff){
	parseLine(buff,-1);
	clauses = new vector<int>[numCs];
	posC= new vector<int>[numVs];
	negC= new vector<int>[numVs];
	for(int i = 0; i < numVs; i++){
		posC[i].push_back(4);
		posC[i].push_back(0);
		posC[i].push_back(0);
		posC[i].push_back(0);
	}
	for(int i = 0; i < numVs; i++){
		negC[i].push_back(4);
		negC[i].push_back(0);
		negC[i].push_back(0);
		negC[i].push_back(0);
	}
	clauseT.reserve(numVs);
}
void parseLine(string line,int indexC){
	char* str = strdup(line.c_str());
    const char s[2] = " ";
   if(inter){
	   if( indexC == -1){
		    /*strtok(str, s);
			numVs = atoi(strtok(NULL, s))+1;
			numV1 = atoi(strtok(NULL, s))+1;
			numCs = atoi(strtok(NULL, s));
			numC1 = atoi(strtok(NULL, s));
			numCc = atoi(strtok(NULL, s));
			return;
			*/
			strtok(str, s);
			numVs = atoi(strtok(NULL, s))+1;
			strtok(NULL, s);
			numV1 = 0;
			numCs = atoi(strtok(NULL, s));
			numC1 = 0;
			numCc = 0;
			return;
	   }
   }// for partition file;*/
   else{
		   if( indexC == -1){
			strtok(str, s);
			strtok(NULL, s);
			numVs = atoi(strtok(NULL, s))+1;
			numV1 = numVs/2;
			numCs = atoi(strtok(NULL, s));
			numC1 = 0;
			numCc = numCs;
			return;
	   }
    }
    int lit;
    int size;
    int p;
    if(indexC < numC1) p = 1;
    else{
    	if(indexC < numCc) p = 2;
    	else p = 3;
    }
    char* token = strtok(str, s);
    while(token != NULL){
		if(*token== '-'){
			lit = atoi(token);
			clauseT.push_back(lit);
		    negC[-lit][p]++;
			token = strtok(NULL, s);
			size++;
			continue;
		}
		if(*token == '0'){
			clauses[indexC] = clauseT;
			clauseT.clear();
		    return;
		}
		lit = atoi(token);
		clauseT.push_back(lit);
	    posC[lit][p]++;
	    size++;
		token = strtok(NULL, s);
    }
	perror("a clause line does not terminates in parseLine");
	exit(EXIT_FAILURE);
}
template<class T>
void Process<T>::printOptions(){
	printf("localSAT options: \n");
	cout<<"c tabu_flag: "<<tabu_flag<<endl;
	cout<<"c maxSteps: "<<maxSteps<<endl;
	cout<<"c seed: "<<seed<<endl;
	cout<<"c fct: "<<fct<<endl;
	cout<<"c cct: "<<cct<<endl;
	cout<<"c cb: "<<cb<<endl;
	cout<<"c eps: "<<eps<<endl;
	cout<<"c lct: "<<lct<<endl;
	switch(fct){
	case 0:{
		cout<<"c polynomial function"<<endl;
		cout<<"c eps: "<<eps<<endl;
		cout<<"c cb: "<<cb<<endl;
		cout<< "pow((eps+break),-cb)" << endl;
		break;
		  }
	default:{
		cout<<"c exponential function"<<endl;
		cout<<"c cb: "<<cb<<endl;
		cout<< "pow(cb,-break)"<< endl;
		break;
		   }
	}
}
void printVariables(){
	cout<< "Variables "<< ": " <<endl ;
   	for(int i = 1; i < numVs; i++){
   		cout<< "Variable "<<i << ": " <<endl ;
   		printVector(posC[i]);
   		printVector(negC[i]);
   	}
}
void printClauses(){
	cout<< "Clauses "<< ": " << endl ;
   	for(int i = 0; i < numCs; i++){
   		cout<< "Clause "<< i<< ": " ;
   		printVector(clauses[i]);
   	}
}
template<class T>
void Process<T>::printAssignment(){
	cout<< "v ";
	for(int i = 1; i < numVs; i++){
		if(assign[i]) cout <<i<<" ";
		else cout << -i<<" ";
	}
	cout <<endl ;
}
template<class T>
void Process<T>::printUnsatCs(){
	cout<< "Unsatisfied clauses ";
	printVector(unsat[0]);
	cout <<endl ;
	printVector(unsat[1]);
	cout <<endl ;
	printVector(unsat[2]);
	cout <<endl ;
}
template<class T>
void Process<T>::printNumP(){
	cout<< "numP: ";
	for(int i = 0; i < numCs; i++){
		cout << numP[i]<< " ";
	}
	cout<<endl;
}




void initialAssignment(){
	for(int i = 0; i < numVs; i++){
		posC[i][1] = posC[i][1]+posC[i][0];
		posC[i][2] =  posC[i][1]+posC[i][2];
		posC[i][3] = posC[i][2]+posC[i][3];
		negC[i][1] = negC[i][1]+negC[i][0];
		negC[i][2] =  negC[i][1]+negC[i][2];
		negC[i][3] = negC[i][2]+negC[i][3];
		posC[i].reserve(posC[i][3]);
		negC[i].reserve(negC[i][3]);
	}
	for(int j = 0; j < numCs; j++){
		for (std::vector<int>::const_iterator i = clauses[j].begin(); i != clauses[j].end(); ++i){
			if(*i < 0) negC[-(*i)].push_back(j);
			else  posC[(*i)].push_back(j);
		}
	}
}
template<class T>
void Process<T>::biasAssignment(){
	for(int i = 0; i < numVs; i++){
			if(posC[i][3] > negC[i][3]){
				assign[i] = true;
			}
			else{
				assign[i] = false;
			}
	}
	setAssignment();
}
template<class T>
void Process<T>::randomBiasAssignment(){
	int sum;
	for(int i = 0; i < numVs; i++){
		sum = posC[i][3] +negC[i][3];
		if(sum == 0){
			assign[i] = true;
		}
		else{
			assign[i] = ((this->*randINT)()%sum)<posC[i][3];
		}
	}
	setAssignment();
}
template<class T>
void Process<T>::randomAssignment(){
   	for(int j = 0; j < numVs; j++){
   		assign[j] = ((this->*randINT)()%2 ==1);
   	}
    setAssignment();
}

template<class T>
void Process<T>::setAssignment(){
   	for(int i = 0; i < numCs; i++){
   		numP[i] = 0;
   	}
	if( tabu_flag && (this->*randINT)()%100<cct){
		for(int i = 0; i < numVs; i++){
			tabuS[i] =0;
		}
	}
   	for(int j = 0; j < numVs; j++){
		if(assign[j] == false){
	   		for (std::vector<int>::const_iterator i = negC[j].begin()+4; i != negC[j].end(); ++i){
	   			numP[*i]++;
	   		}
		}
		else{
			for (std::vector<int>::const_iterator i = posC[j].begin()+4; i != posC[j].end(); ++i){
	   			numP[*i]++;
			}
   		}
   	}
   	for(int i = 0; i < numC1; i++){
   		if(numP[i] == 0){
   			numUnsat++;
   			unsat[0].push_back(i);
   		}
   	}
   	for(int i = numC1; i < numCc; i++){
   		if(numP[i] == 0){
 			numUnsat++;
   			unsat[1].push_back(i);
   		}
   	}
   	for(int i = numCc; i < numCs; i++){
   		if(numP[i] == 0){
 			numUnsat++;
   			unsat[2].push_back(i);
   		}
   	}
}
template<class T>
void Process<T>::solvePart(int index){
	assert(numUnsat == (unsat[0].size()+unsat[1].size()+unsat[2].size()));
	vector<int>& unsatCs = unsat[index];
	for(unsigned int j = 0; j < maxSteps; j++){
		if (unsatCs.size()== 0){
			return;
		}
		int size = unsatCs.size();
		int randC = (this->*randINT)()%size;
		int flipCindex = unsatCs[randC];
		while(numP[flipCindex] > 0){
			unsatCs[randC]=unsatCs.back();
			unsatCs.pop_back();
			size--;
 			numUnsat--;
 			assert(numUnsat == (unsat[0].size()+unsat[1].size()+unsat[2].size()));
			if(size == 0) return;
			randC = (this->*randINT)()%size;
			flipCindex = unsatCs[randC];
		}
		int flipLindex = getFlipLiteral(flipCindex);
		unsatCs[randC]=unsatCs.back();
		unsatCs.pop_back();
		numUnsat--;
		flip(flipLindex);
		if(tabu_flag) tabuS[abs(flipLindex)]++;
	}

	return;
}
template<class T>
void Process<T>::optimal(){
	int rct;
	//check(get first unsat part)
	//solveP0
	//solveP1
	//solveC
	while(true){
			if (numUnsat == 0){
				assert(numUnsat == (unsat[0].size()+unsat[1].size()+unsat[2].size()));
			/*	#pragma omp critical
				{
				test();
				}
			*/
				cout<< "s SATISFIABLE"<< endl;
				//printAssignment();
				//abort();
				return;
			}
			solvePart(0);
			assert(numUnsat == (unsat[0].size()+unsat[1].size()+unsat[2].size()));
			solvePart(2);
			assert(numUnsat == (unsat[0].size()+unsat[1].size()+unsat[2].size()));
			solvePart(1);
			assert(numUnsat == (unsat[0].size()+unsat[1].size()+unsat[2].size()));
	}
}

template<class T>
int Process<T>::getFlipLiteral(int cIndex){
	vector<int>&  vList = clauses[cIndex];
	int j=0,bre,min= numCs+1;
	double sum=0,randD;
	int greedyLiteral = 0, randomLiteral;
	for (std::vector<int>::const_iterator i = vList.begin(); i != vList.end(); ++i){
		bre = computeBreakScore(*i);
		if(bre == 0){
			if (numUnsat < ((this->*randINT)()%numCs)){
				return *i;
			}
		}
		if(tabu_flag &&bre == 0 && tabuS[abs(*i)] == 0) return *i;
		if(bre < min){
			min = bre;
			greedyLiteral = *i;
		}
		if(bre < numCs){
		sum+= lookUpTable[bre];
		}
		else{
		sum+=(this->*Process::lookUp)(bre);
		}
		probs[j]= sum;
		j++;
	}
	randD = ((double)(this->*randINT)()/RAND_MAX)*sum;
	assert(randD >= 0);
	for(int i = 0; i < j;i++){
		if(probs[i]< randD){
			continue;
		}
		randomLiteral= vList[i];
		break;
	}
	if(tabu_flag &&tabuS[abs(greedyLiteral)] < tabuS[abs(randomLiteral)]){
		return greedyLiteral;
	}
	return randomLiteral;
}
template<class T>
void Process<T>::flip(int literal){
	std::vector<int>::const_iterator i;
	if(literal > 0){
   		for (i = negC[literal].begin()+4; i != negC[literal].end(); ++i){
   			numP[*i]--;
   			if(numP[*i] == 0) push(*i);
   		}
		for (i = posC[literal].begin()+4; i != posC[literal].end(); ++i){
   			numP[*i]++;
		}

		assign[literal] = true;
	}
	else{
   		for (i = negC[-literal].begin()+4; i != negC[-literal].end(); ++i){
   			numP[*i]++;
   		}
		for (i = posC[-literal].begin()+4; i != posC[-literal].end(); ++i){
   			numP[*i]--;
   			if(numP[*i] == 0) push(*i);
		}
		assign[-literal]= false;
	}
}
template<class T>
void Process<T>::test(){
	ifstream fp;
	fp.open(fileName,std::ios::in);
	if(!fp.is_open()){
		perror("read file fails");
		exit(EXIT_FAILURE);
	}
	string buff;
	char head;
   	getline(fp,buff);
   	while(!fp.eof()){
   		if(buff.empty()) break;
		head =buff.at(0);
		if(inter){
			if(head == 'c'){
				break;
			}
		}
		else{
			if(head == 'p'){
				break;
			}
		}
	  getline(fp,buff);
	}
   	while(!fp.eof()){
		getline(fp,buff);
		if(buff.empty()) break;
		testLine(buff);
   	}
   	cout<< "tested" << endl;
}
template<class T>
void Process<T>::testLine(string line){
	char* str = strdup(line.c_str());
    const char s[2] = " ";
    int lit;
    int numT=0;
    char* token = strtok(str, s);
    while(token != NULL){
		if(*token== '-'){
			lit = atoi(token);
			if(assign[-lit] == false) numT++;
			token = strtok(NULL, s);
			continue;
		}
		if(*token == '0'){
			if(numT == 0){
				perror("TEST FAILURE");
				exit(EXIT_FAILURE);
			}
		    return;
		}
		lit = atoi(token);
		if(assign[lit] == true) numT++;
		token = strtok(NULL, s);
    }
    cout<< line;
	perror("a clause line does not terminates in testLine");
	exit(EXIT_FAILURE);

}
template<class T>
int Process<T>::computeBreakScore(int literal){
    int score = 0;
    int aIndex = abs(literal);
    vector<int>& occList =(literal < 0)? posC[aIndex] :negC[aIndex];
    for(std::vector<int>::const_iterator i = occList.begin()+4; i != occList.end(); ++i) {
        if (numP[*i]== 1) {
            score++;
        }
    }
	//cout<< "out make "<<endl;
    return score;
}
template<class T>
double Process<T>::func_exp(int literal){
	return pow(cb,-computeBreakScore(literal));
}
template<class T>
double Process<T>::func_poly(int literal){
	return pow((eps+computeBreakScore(literal)),-cb);
}


void printVector(vector<int>& vec){
	for (std::vector<int>::const_iterator i = vec.begin(); i != vec.end(); ++i){
		cout << *i << ' ';
	}
	cout << endl;
}

/*print the manual of the program*/
void printUsage(){
	printf("\n");
	printf("\n");
	printf("\n");
	printf("\n--------localSAT, a stochastic local search solving SAT problem-----------------\n");
	printf("\nautor: Guangping Li\n");
	printf("references: Engineering a lightweight and efficient local search sat solver \n");

	printf("\n-------------------Usage--------------------------------------------------------\n");
	printf("\n");
	printf("./locolSAT  <DIMACS CNF instance> [options]\n");
	printf("\n--------------------------------------------------------------------------------\n");

	printf("\n");
	printf("\n-------------------Options------------------------------------------------------\n");
	printf("\n");
	printf("**for the potential flipping function:\n");
	printf("   --function, -f: 0 =  polynomial; 1 = exponential with make; 2 = exponential only with break: 3 = equal  [default = 0]\n");
	printf("   --eps, -e <double value> [0.0,1.0] :  [default = 1.0]\n");
	printf("   --cb, -b <double value>  [2.0,10.0] : constant for break\n");
	printf("   --cm, -m <double value>  [-2.0,2.0] : constant for make\n");
	printf("   --w, -w <double value>  [0.0,1.0] : [default = 0.3] the probability used only in the WALKSAT\n");
	printf("\n");
	printf("**for the process:\n");
	printf("   --caching, -c  : 0 =  no caching ; 1 = caching  for scores\n");
	printf("   --rounds, -r <long long int_value>,   : maximum number of tries \n");
	printf("   --flips, -p  <long long int_value>,   : maximum number of flips per try \n");
	printf("   --seed, -s : seed used for Randomness\n");
	printf("\n");
	printf("**further options:\n");
	printf("   --output, -o : output the final assignment\n");
	printf("   --help, -h : output this help\n");
	printf("\n");
	printf("---------------------------------------------------------------------------------\n");
}

template<class T>
void Process<T>::initLookUpTable_exp(){

	lookUpTable = (double*)malloc(sizeof(double) * numCs);
	for(int i = 0; i < numCs;i++){
		lookUpTable[i] = pow(cb,-i);
	}
}
template<class T>
void Process<T>::initLookUpTable_poly(){
	lookUpTable = (double*)malloc(sizeof(double) * numCs);
	for(int i = 0; i < numCs;i++){
		lookUpTable[i] = pow((eps+i),-cb);
	}
}
template<class T>
double Process<T>::LookUpTable_exp(int bre){
	return pow(cb,-bre);
};
template<class T>
double Process<T>::LookUpTable_poly(int bre){
	return pow((eps+bre),-cb);
};
template<class T>
int Process<T>::randI(){
	return distribution(generator);
};

template<class T>
int Process<T>::randI2(){
	return rand();
};
template<class T>
void Process<T>::push(int cIndex){
	assert(numUnsat == (unsat[0].size()+unsat[1].size()+unsat[2].size()));
	numUnsat++;
	if(cIndex < numC1) unsat[0].push_back(cIndex);
	else{
		if(cIndex < numCc) unsat[1].push_back(cIndex);
		else unsat[2].push_back(cIndex);

	}
	assert(numUnsat == (unsat[0].size()+unsat[1].size()+unsat[2].size()));
}


