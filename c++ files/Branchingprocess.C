#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TFormula.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>
#include <TF1.h>
#include <string>
#include <TObject.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>



#endif

using namespace std;
void branchingprocess(double p, int & active, int & generations);

int branching(int & active, int &generations, double p, vector <int> & gens);

void sim(int iterations, double p){
	
	
	ofstream fout("outputbranching.csv");

	for (int i = 0; i < iterations; i ++){

		int active = 1;
		int generations = 0;
	
		branchingprocess(p,active, generations);

		fout << active << ",";
		fout << generations << "\n";

	}
fout.close();
}


int branching(int & active, int & generations, double p, vector <int> & gens){
	bool debug = kFALSE;
	int branches = 2;
	int maxgenerations = 100000;
	for (int i = 0; i < branches; i ++){
		
		int gen = generations;
		double n = gRandom->Uniform(0,1);
		if (debug) cout << "Evaluating node number "<< i <<  endl;
		if (debug) cout << "Generations = " << gen << endl;

		if (n < p && gen < maxgenerations){
			
			gen ++;
			active = active + 2;
			if (debug) printf("%d nodes are created \n", branches);
			if (debug) cout <<"Generation number " << gen <<  " is created" << endl;
			

			branching(active,gen,p,gens);
	
		}

		else {
			if (debug) cout <<"This node has died at generation number " << gen<< endl;
			gens.push_back(gen);
			

		}
		

	}


	if (debug) cout << "Max generation is " << *max_element(gens.begin(), gens.end()) << endl;
	return  *max_element(gens.begin(), gens.end());
}


void branchingprocess(double p, int & active, int & generations){
	bool debug = kFALSE;

	if (debug) cout <<"First generation" << endl;
	int num;
	int branches = 2;
	vector< int> gens;
	if (debug) cout << "Generations = " << generations << endl;
	double n = gRandom->Uniform(0,1);
	if (n < p){

		generations++;
		active = active + branches;

		if (debug) printf("%d nodes are created \n", branches);
		if (debug) cout << "Generation number " << generations <<  " is created" << endl;
		if (debug) cout << "Generations = " << generations << endl;
		

		num = branching(active,generations,p,gens);
		generations = num;

	
	}

	else{

	if (debug) cout << "Avalanche has ended immediately. Generations = " << generations << endl;
	}
	
}

