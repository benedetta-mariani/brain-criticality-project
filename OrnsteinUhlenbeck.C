#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TFormula.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TF1.h>
#include <string>
#include <TObject.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>

#endif

using namespace std;


double* DiffusionCoefficient(double dt, double T, double tau, double D);
double* OrnsteinUhlebeck(double dt, double T, double tau, double* D);


void Main(int nunits, int ntime){ 

	double dt = 0.001;
	double T = ntime;
	double tau = 10;
	double Di = 1;
	int N = int(T/dt);

	double tau2 = 0.08;

	ofstream fout("output.csv");
	double* diff;
	diff = DiffusionCoefficient(dt, T, tau, Di);
	
	for (int i = 0; i < nunits; i ++){
		
		double* poi;
		poi = OrnsteinUhlebeck(dt, T, tau2, diff);
		fout << *(poi);
		for (int s = 1; s < N; s ++){
			fout << ", "<< *(poi+s);
		}

		delete []poi;
		fout << " \n ";

	}

	delete []diff;
	fout.close();
}



double* OrnsteinUhlebeck(double dt, double T, double tau, double* D) {
	
	int N = int(T/dt);
	double *x = new double[N];
	x[0] = 0;
	for (int t = 0; t < N; t ++){
		x[t+1] = x[t] -x[t]*dt/tau + TMath::Sqrt(*(D+t)*dt)*gRandom->Gaus(0,1);
	}

	return x;
	
	
}

double* DiffusionCoefficient(double dt, double T, double tau, double D) {

	int N = int(T/dt);
	double *x = new double[N];	
	x[0] = 0;
	for (int t = 0; t < N; t ++){
		x[t+1] = x[t] -x[t]*dt/tau + TMath::Sqrt(D*dt)*gRandom->Gaus(0,1);
	}

	for (int i = 0; i < N; i ++){
		if (x[i] < 0.2){
			x[i] = 0.2;
		}
	}
	return  x;
	
	
}


    
	


	




