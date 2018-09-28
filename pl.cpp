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
#pragma omp parallel num_threads(2)
 {
	const vector<int> setI =setII[0];
	const vector<double>& setD = setDD[0];
	switch(omp_get_thread_num()){
	case 0:{
		Process<minstd_rand0> process;
		process.optimal();
		break;
	}
	case 1:{
		Process<minstd_rand> process;
		process.optimal();
		break;
	}
	case 2:{
		Process<mt19937> process;
		process.optimal();

		break;
	}
	case 3:{
		Process<mt19937_64> process;
		process.optimal();
		break;
	}
	case 4:{
		Process<ranlux24_base> process;
		process.optimal();
		break;
	}
	case 5:{
		Process<ranlux48_base> process;
		process.optimal();
		break;
	}
	case 6:{
		Process<ranlux24> process;
		process.optimal();
		break;
	}
	case 7:{
		Process<ranlux48> process;
		process.optimal();
		break;
	}
	case 8:{
		Process<knuth_b> process;
		process.optimal();
		break;
	}
	case 9:{
		Process<default_random_engine> process;
		process.optimal();
		break;
	}
	case 10:{
		Process<minstd_rand0> process;
		process.optimal();
		break;
	}
}
}
 testPart(0,assignG);
 testPart(1,assignG);
 testPart(2,assignG);
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
Process<T>::Process():distribution(0, INT_MAX){
	//set the parameters
	   // set tabuS
	clauseQ.reserve(maxL);
	if(omp_get_thread_num() == 0){
		srand(seed);
		randINT = &Process::randI2;
	}
	else{
		randINT = &Process::randI;
		generator.seed(seed);

	}
		tabuS = (int*) malloc(sizeof(int) * numVs);
		for(int i = 0; i < numVs; i++){
			tabuS[i] = 0;
		}
	numP = (int*) malloc(sizeof(int) * numCs);
	probs = (double*)malloc(sizeof(double) * numVs);
	assign = (bool*)malloc(sizeof(bool) * numVs);
	unsat.reserve(numCs);

	if(maxL <= 3){
		initAssignment = &Process::randomAssignment;
		cb = 2.06;
		eps = 0.9;
		fct = 0;
		setAssignment =&Process::setAssignment3;
		getFlipLiteral =&Process::getFlipLiteral3;
		flip = &Process::flip3;
	}
	else if (maxL <=4){
		initAssignment = &Process::randomAssignment;
		cb = 2.85;
		setAssignment =&Process::setAssignment3;
		getFlipLiteral =&Process::getFlipLiteral3;
		flip = &Process::flip3;
	}
	else if(maxL <=5){
		initAssignment = &Process::biasAssignment;
		cb = 3.7;
		breaks = (int*) malloc(sizeof(int) * numVs);
		critVar = (int*) malloc(sizeof(int) * numCs);
		setAssignment =&Process::setAssignment57;
		getFlipLiteral =&Process::getFlipLiteral57;
		flip = &Process::flip57;
	}
	else if(maxL <= 6){
		breaks = (int*) malloc(sizeof(int) * numVs);
		critVar = (int*) malloc(sizeof(int) * numCs);
		initAssignment = &Process::biasAssignment;
		cb = 5.1;
		setAssignment =&Process::setAssignment57;
		getFlipLiteral =&Process::getFlipLiteral57;
		flip = &Process::flip57;
	}
	else{
		breaks = (int*) malloc(sizeof(int) * numVs);
		critVar = (int*) malloc(sizeof(int) * numCs);
		initAssignment = &Process::biasAssignment;
		cb = 5.4;
		setAssignment =&Process::setAssignment57;
		getFlipLiteral =&Process::getFlipLiteral57;
		flip = &Process::flip57;
	}
	switch (fct){
	case 0:initLookUpTable_poly();
			lookUp =&Process::LookUpTable_poly;
			break;
	default:initLookUpTable_exp();
			lookUp =&Process::LookUpTable_exp;
			break;
	}
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
	assignG = (bool*)malloc(sizeof(bool) * numVs);
	static const int arr[] = {4,0,0,0};
	for(int i = 0; i < numVs; i++){
		posC[i]= {4,0,0,0};
		/*posC[i].push_back(4);
		posC[i].push_back(0);
		posC[i].push_back(0);
		posC[i].push_back(0);*/
	}
	for(int i = 0; i < numVs; i++){
		negC[i]= {4,0,0,0};
		/*negC[i].push_back(4);
		negC[i].push_back(0);
		negC[i].push_back(0);
		negC[i].push_back(0);*/
	}
	clauseT.reserve(numVs);
}
void parseLine(string line,int indexC){
	char* str = strdup(line.c_str());
    const char s[2] = " ";
   if(inter){
	   if( indexC == -1){
		    strtok(str, s);
			numVs = atoi(strtok(NULL, s))+1;
			numV1 = atoi(strtok(NULL, s))+1;
			numCs = atoi(strtok(NULL, s));
			numC1 = atoi(strtok(NULL, s));
			numCc = atoi(strtok(NULL, s));
			return;
	/*
			strtok(str, s);
			numVs = atoi(strtok(NULL, s))+1;
			strtok(NULL, s);
			numV1 = numVs;
			numCs = atoi(strtok(NULL, s));
			numC1 = 0;
			numCc = numCs;
			return;
	*/
	   }
   }// for partition file;*/
   else{
		   if( indexC == -1){
			strtok(str, s);
			strtok(NULL, s);
			numVs = atoi(strtok(NULL, s))+1;
			numV1 = 0;
			numCs = atoi(strtok(NULL, s));
			numC1 = 0;
			numCc = 0;
			return;
	   }
    }
    int lit;
    int size = 0;
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
			if(size > maxL) maxL = size;
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
	printVector(unsat);

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
void Process<T>::biasAssignment(int partition){
	int vs,ve, start,end;
	switch(partition){
	case 0: vs = 0; ve = numV1; start = 0; end = 1; break;
	case 1: vs = numV1; ve = numVs; start = 2; end = 3;break;
	default:vs = 0; ve = numVs; start = 0; end = 3;
	}
	for(int i = vs; i < ve; i++){
			if((posC[i][end]-posC[i][start]) > (negC[i][end]-negC[i][start])){
				assign[i] = true;
			}
			else{
				assign[i] = false;
			}
	}
	(this->*setAssignment)(partition);
}
template<class T>
void Process<T>::randomAssignment(int partition){
	int vs,ve;
	switch(partition){
	case 0: vs = 0; ve = numV1; break;
	case 1: vs = numV1; ve = numVs; break;
	default:vs = 0; ve = numVs;
	}
   	for(int j = vs; j < ve; j++){
   		assign[j] = ((this->*randINT)()%2 ==1);
   	}
	(this->*setAssignment)(partition);
}

template<class T>
void Process<T>::setAssignment3(int partition){
	unsat.clear();
	int cs, ce;
	switch(partition){
	case 0:cs = 0; ce = numC1;break;
	case 2:cs = numCc; ce = numCs;break;
	default: cs = 0; ce = numCs;
	}
   	for(int i = cs; i < ce; i++){
   		numP[i] = 0;
   	}
   	for(int j = cs; j < ce; j++){
   		for (std::vector<int>::const_iterator i = clauses[j].begin(); i != clauses[j].end(); ++i){
   			if(((*i) > 0 && (assign[*i]== true))||((*i) < 0 && assign[-*i]== false)){
   				numP[j]++;
   			}
   		}
   		if(numP[j] == 0){
   			unsat.push_back(j);
   		}
   	}

}
template<class T>
void Process<T>::setAssignment57(int partition){
	unsat.clear();
	int cs, ce, vs, ve;
	switch(partition){
	case 0:cs = 0; ce = numC1;vs = 0; ve = numV1;break;
	case 2:cs = numCc; ce = numCs;vs = numV1; ve = numVs;break;
	default: cs = 0; ce = numCs;vs = 0; ve = numVs;
	}
   	for(int i = cs; i < ce; i++){
   		numP[i] = 0;
   	}
	for (int i = vs; i <= ve; i++) {
		breaks[i] = 0;
	}
	int critV = 0;
   	for(int j = cs; j < ce; j++){
   		for (std::vector<int>::const_iterator i = clauses[j].begin(); i !=clauses[j].end(); ++i){
   			if(((*i) > 0 && assign[*i]== true)||((*i) < 0 && assign[-*i]== false)){
   				numP[j]++;
   				critV = *i;
   			}
   		}
   		if(numP[j] == 0){
   			unsat.push_back(j);
   		}
   		else{
   			if(numP[j] == 1){
   				critV = abs(critV);
   				critVar[j] = critV;
   				breaks[critV]++;
   			}
   		}
   	}
}

template<class T>
void Process<T>::solve(){
	while(true){
		if(sat) return;
		if (!sat && unsat.size()== 0){
			#pragma omp critical
			{
			for(int i = 0; i < numVs; i++){
						assignG[i] = assign[i];
			}
			}
			sat = true;
			return;
		}
		int size = unsat.size();
		int randC = (this->*randINT)()%size;
		int flipCindex = unsat[randC];
		while(numP[flipCindex] > 0){
			unsat[randC]=unsat.back();
			unsat.pop_back();
			size--;
			if(sat) return;
			if (!sat && size == 0){
				sat = true;
				#pragma omp critical
				{
					for(int i = 0; i < numVs; i++){
						assignG[i] = assign[i];
					}
				}
				return;
				}
			randC = (this->*randINT)()%size;
			flipCindex = unsat[randC];
		}
		if(sat) return;
		int flipLindex = (this->*getFlipLiteral)(flipCindex,-1);
		unsat[randC]=unsat.back();
		unsat.pop_back();
		if(sat) return;
		(this->*flip)(flipLindex,-1);
		tabuS[abs(flipLindex)]++;
	}
}
template<class T>
void Process<T>::solvePart(int index){
	bool& sat = (index == 0)? sat0: sat2;
	while(true){
		if(sat) return;
		if (!sat && unsat.size()== 0){
			if(index == 0){
				#pragma omp critical
				{
					for(int i = 0; i < numV1; i++){
						assignG[i] = assign[i];
					}
				}
			}
			else{
				#pragma omp critical
				{
					for(int i = numV1; i < numVs; i++){
						assignG[i] = assign[i];
					}
				}
			}
			sat = true;
			return;
		}
		int size = unsat.size();
		int randC = (this->*randINT)()%size;
		int flipCindex = unsat[randC];
		while(numP[flipCindex] > 0){
			unsat[randC]=unsat.back();
			unsat.pop_back();
			size--;
			if(sat) return;
			if (!sat && size == 0){
				if(index == 0){
					#pragma omp critical
					{
						for(int i = 0; i < numV1; i++){
							assignG[i] = assign[i];
						}
					}
				}
				else{
					#pragma omp critical
					{
						for(int i = numV1; i < numVs; i++){
							assignG[i] = assign[i];
						}
					}
				}
				sat = true;
				return;
			}
			randC = (this->*randINT)()%size;
			flipCindex = unsat[randC];
		}
		if(sat) return;
		int flipLindex = (this->*getFlipLiteral)(flipCindex,index);
		unsat[randC]=unsat.back();
		unsat.pop_back();
		if(sat) return;
		(this->*flip)(flipLindex,index);
		tabuS[abs(flipLindex)]++;
	}
}
template<class T>
void Process<T>::optimal(){
	bool odd = (omp_get_thread_num()%2 == 0);
	bool& satF = odd? sat0:sat2;
	bool& satS = odd? sat2:sat0;
	int First = odd? 0:2;
	int Second = odd?2:0;
	if(!satF){
		(this->*initAssignment)(First);
		solvePart(First);
		if(unsat.size() != 0){
			if(First == 0){
				for(int i =0; i < numV1; i++){
					assign[i] = assignG[i];
				}
			}
			else{
				for(int i =numV1; i < numVs; i++){
					assign[i] = assignG[i];
				}

			}
		}
	}
	else{
		if(First == 0){
			for(int i =0; i < numV1; i++){
				assign[i] = assignG[i];
			}
		}
		else{
			for(int i =numV1; i < numVs; i++){
				assign[i] = assignG[i];
			}

		}

	}
	if(!satS){
		(this->*initAssignment)(Second);
		solvePart(Second);
		if(unsat.size() != 0){
			if(Second == 0){
				for(int i =0; i < numV1; i++){
					assign[i] = assignG[i];
				}
			}
			else{
				for(int i =numV1; i < numVs; i++){
					assign[i] = assignG[i];
				}

			}

		}

	}
	else{
		if(Second == 0){
			for(int i =0; i < numV1; i++){
				assign[i] = assignG[i];
			}
		}
		else{
			for(int i =numV1; i < numVs; i++){
				assign[i] = assignG[i];
			}

		}

	}
	(this->*setAssignment)(-1);
	solve();
}

template<class T>
int Process<T>::getFlipLiteral3(int cIndex, int partition){
	clauseQ.clear();
	vector<int>&  vList = clauses[cIndex];
	int j=0,bre,min= numCs+1;
	double sum=0,randD;
	int greedyLiteral = 0, randomLiteral;
	switch(partition){
	case 0: computeBreak = &Process::computeBreakScore0; break;
	case 2: computeBreak = &Process::computeBreakScore2; break;
	case -1:computeBreak = &Process::computeBreakScore; break;
	}
	for (std::vector<int>::const_iterator i = vList.begin(); i != vList.end(); ++i){
		bre = (this->*Process::computeBreak)(*i);
		if(bre == 0){
			clauseQ.push_back(*i);
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
	int cS = clauseQ.size();
	if(cS > 0){
		double temp = noise * (double)flipCount/numVs;
		int index = (this->*randINT)()%cS;
		for(int i =0; i < cS; i++){
			greedyLiteral = clauseQ[index];
			if(tabuS[abs(greedyLiteral)] < temp){
				return greedyLiteral;
			}
			index = (index+1)%cS;
		}

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
	return randomLiteral;
}
template<class T>
int Process<T>::getFlipLiteral57(int cIndex, int partition){
	clauseQ.clear();
	vector<int>&  vList = clauses[cIndex];
	int j=0,bre,min= numCs+1;
	double sum=0,randD;
	int greedyLiteral = 0, randomLiteral;
	switch(partition){
	case 0: computeBreak = &Process::computeBreakScore0; break;
	case 2: computeBreak = &Process::computeBreakScore2; break;
	case -1:computeBreak = &Process::computeBreakScore; break;
	}
	for (std::vector<int>::const_iterator i = vList.begin(); i != vList.end(); ++i){
//		bre = breaks[abs(*i)];
		bre = (this->*Process::computeBreak)(*i);
		if(bre == 0){
			clauseQ.push_back(*i);
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
	int cS = clauseQ.size();
	if(cS > 0){
		double temp = noise * (double)flipCount/numVs;
		int index = (this->*randINT)()%cS;
		if(cS< numCs){
			temp= temp*lookUpTable[cS];
		}
		else{
		temp=  temp*(this->*Process::lookUp)(cS);
		}
		for(int i =0; i < cS; i++){
			greedyLiteral = clauseQ[index];
			if(tabuS[abs(greedyLiteral)] < temp){
				return greedyLiteral;
			}
			index = (index+1)%cS;
		}

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
	return randomLiteral;
}
template<class T>
void Process<T>::flip3(int literal,int partition){
	flipCount++;
    int aIndex = abs(literal);
    vector<int>& occList =(literal < 0)? posC[aIndex] :negC[aIndex];
    vector<int>& deList =(literal < 0)? negC[aIndex] : posC[aIndex] ;
    int start, end;
    int startD, endD;
    switch(partition){
    case -1:start = occList[0];end = occList[3]; startD = deList[0];endD = deList[3]; break;
    case 0:start = occList[0]; end = occList[1];startD = deList[0]; endD = deList[1];break;
    case 2:start = occList[2]; end = occList[3];startD = deList[2]; endD = deList[3];break;
    }
	for (int i = start; i <end; ++i){
		numP[occList[i]]--;
		if(numP[occList[i]] == 0){
	 	   unsat.push_back(occList[i]);
		}
	}
	for (int i = startD; i <endD; ++i){
		numP[deList[i]]++;
	}

	if(literal > 0)assign[literal] = true;
	else assign[-literal] = false;
}
template<class T>
void Process<T>::flip57(int literal,int partition){
	flipCount++;
    int aIndex = abs(literal);
    assign[aIndex] = (literal > 0);
    vector<int>& occList =(literal < 0)? posC[aIndex] :negC[aIndex];
    vector<int>& deList =(literal < 0)? negC[aIndex] : posC[aIndex] ;
    int start, end;
    int startD, endD;
    int cla; int num;
    switch(partition){
    case -1:start = occList[0];end = occList[3]; startD = deList[0];endD = deList[3]; break;
    case 0:start = occList[0]; end = occList[1];startD = deList[0]; endD = deList[1];break;
    case 2:start = occList[2]; end = occList[3];startD = deList[2]; endD = deList[3];break;
    }
	for (int i = start; i <end; ++i){
		cla = occList[i];
		num = numP[cla];
		if(num == 1){
			breaks[aIndex]--;
			unsat.push_back(cla);
		}
		else if(num == 2){
			int aj;
			for(std::vector<int>::const_iterator j = clauses[cla].begin(); j != clauses[cla].end(); ++j){
				aj = abs(*j);
				if (((*j)==aj)== assign[aj]) {
					critVar[cla] = aj;
					breaks[aj]++;
					break;
				}
			}

		}
		numP[cla]--;
	}
	for (int i = startD; i <endD; ++i){
		cla = deList[i];
		num = numP[cla];
		if(num== 0){
			critVar[cla] = aIndex;
			breaks[aIndex]++;
		}
		else if(num == 1){
			breaks[critVar[cla]]--;
		}
		numP[cla]++;
	}
}
void testPart(int partition, bool* assignG){
	int num;
	switch(partition){
		case 0:num = numC1;break;
		case 1:num = numCc-numC1;break;
		case 2:num = numCs - numCc; break;
		assert(false);
	}
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
   	int line = 0;
   	if(partition == 0){
		while(!fp.eof()){
			if(line == num){
				cout<< partition << " tested" << endl;
				return;
			}
			getline(fp,buff);
			if(buff.empty()) continue;
			testLine(buff);
			line++;
		}
   	}
   	if(partition == 1){
		while(!fp.eof()){
			if(line == numC1) break;
			getline(fp,buff);
			line++;
   		}
	   	line = 0;
		while(!fp.eof()){
			if(line == num) {
				cout<< partition << " tested" << endl;
				return;
			}
			getline(fp,buff);
			if(buff.empty()) continue;
			testLine(buff);
			line++;
		}
	   	cout<< num <<" 0 line number "<<  line <<endl;
	   	assert(false);
   	}
   	if(partition == 2){
   	   	line = 0;
		while(!fp.eof()){
			if(line == numCc) break;
			getline(fp,buff);
			if(buff.empty()) continue;
			line++;
		}
	   	line = 0;
		while(!fp.eof()){
			if(line == num) {
				cout<< partition << " tested" << endl;
				return;
			}
			getline(fp,buff);
			if(buff.empty()){
				cout<< " "<<line <<"empty line";
				break;
			}
			testLine(buff);
			line++;
		}
	   	cout<<"2 line number "<<  line <<endl;
	   	assert(false);

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
void testLine(string line){
	char* str = strdup(line.c_str());
    const char s[2] = " ";
    int lit;
    int numT=0;
    char* token = strtok(str, s);
    while(token != NULL){
		if(*token== '-'){
			lit = atoi(token);
			if(assignG[-lit] == false) numT++;
			token = strtok(NULL, s);
			continue;
		}
		if(*token == '0'){
			if(numT == 0){
				cout<< line<< endl;
				perror("TEST FAILURE");
				exit(EXIT_FAILURE);
			}
		    return;
		}
		lit = atoi(token);
		if(assignG[lit] == true) numT++;
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
    int start, end;
    start = occList[0];end = occList[3];
    for(int i = start; i < end; i++) {
        if (numP[occList[i]]== 1) {
            score++;
        }
    }
	//cout<< "out make "<<endl;
    return score;
}
template<class T>
int Process<T>::computeBreakScore0(int literal){
    int score = 0;
    int aIndex = abs(literal);
    vector<int>& occList =(literal < 0)? posC[aIndex] :negC[aIndex];
    int start, end;
    start = occList[0]; end = occList[1];
    for(int i = start; i < end; i++) {
        if (numP[occList[i]]== 1) {
            score++;
        }
    }
	//cout<< "out make "<<endl;
    return score;
}
template<class T>
int Process<T>::computeBreakScore2(int literal){
    int score = 0;
    int aIndex = abs(literal);
    vector<int>& occList =(literal < 0)? posC[aIndex] :negC[aIndex];
    int start, end;
    start = occList[2]; end = occList[3];
    for(int i = start; i < end; i++) {
        if (numP[occList[i]]== 1) {
            score++;
        }
    }
	//cout<< "out make "<<endl;
    return score;
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


