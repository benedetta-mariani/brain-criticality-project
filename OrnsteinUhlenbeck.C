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
#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>


#endif

using namespace std;


vector <double> OUGillespieDiffusion(double dt, double T, double tau, double D, double Dstar);
vector<vector<double>> units(int nunits, vector <double>a[], vector<double> diff,double tauD, vector <double> tau, double dt, bool connected, int N);


void Main(int nunits, int T){ 
	double dt = 0.001;
	double tauD = 4;
	double Di = 1;
	vector <double> tau(nunits, 0.);
	for (int i = 0; i < nunits; i ++){
		tau[i] = 0.06;
	} 
	double Dstar = 0.2;
	bool connected = kTRUE;
	int N = int(T/dt);
	int n = int(tauD/dt*100) + N; 

	vector <double> diff;
	vector<double>a[nunits];

	// csv reader
	// Reading the adjacency matrix "a" from a csv file W.csv" 
  	fstream file;

  	file.open("W.csv");
  	string line;
  	int l = 0;
  	while (getline( file, line, '\n')){
	  istringstream templine(line); 
	  string data;
	  while (getline( templine, data, ',')) 
	  {
	    a[l].push_back(atof(data.c_str())); 
	  }
	  l ++;
	}
  	file.close();
  	////////////////////////////////////////////////////////////

  	
  			
	diff = OUGillespieDiffusion(dt, T, tauD, Di, Dstar);
	vector<vector<double>> neuro(nunits, vector<double> (N,0.));
	neuro = units(nunits, a, diff,tauD,tau, dt, connected, N);
	

	// csv writer
	// Writing the resulting signals on a csv file
	ofstream fout("output.csv");
	for (int i = 0; i < nunits; i ++){
		fout << neuro[i][0];
		for (int j = 1; j < N; j ++){
			fout << "," << neuro[i][j];
		}
		fout << "\n";
	}
	fout.close();
	
	///////////////////////////////////////////////	
}




vector <double> OUGillespieDiffusion(double dt, double T, double tau, double D, double Dstar){
	int N = int(T/dt);
	int n = int(tau/dt*100) + N; // to avoid initial transient
	vector <double> x(n,0.);
    	double mu = TMath::Exp(-dt/tau);
    	double sigma = TMath::Sqrt(D*tau/2*(1-pow(mu,2)));
    	x[0] = 0.;
    	for (int t = 0; t < n-1; t ++){
        	x[t+1] = x[t]*mu + sigma*gRandom->Gaus(0,1);  
     	}
	for (int i = 0; i < n; i ++){
		if (x[i] < Dstar){
			x[i] = Dstar;
		}
	}	
	return x;
}




vector<vector<double>> units(int nunits,  vector<double> a[],vector<double> diff, double tauD, vector <double> tau, double dt, bool connected, int N){
	vector<vector<double>> x(nunits, vector<double> (N, 0.));

	double x0;
	for (int i = 0; i < nunits; i ++){
		x0 = 0.;
		x[i][0] = x0;
	}
	
	double number;
	for (int t = 0; t < N-1; t ++){

		for (int i = 0; i < nunits; i ++){
			double value = 0;
			if (connected){
				for (int j = 0; j < nunits; j ++){

					value += a[i][j]*x[j][t]*dt;	
				}
			}
			x[i][t+1] = x[i][t] - x[i][t]*dt/tau[i] + value + TMath::Sqrt(diff[int(tauD*100/dt)+ t]*dt)*gRandom->Gaus(0,1);
			
		}
	}

	return x;

}





