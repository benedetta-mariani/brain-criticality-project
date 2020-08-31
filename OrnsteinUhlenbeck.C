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


double* OrnsteinUhlebeck(double*x,double dt, double T, double tau, double D);
double* OrnsteinUhlebeck2(double*x,double dt, double T, double tau, double* D);




void Main(int nunits, int ntime){ 

	double dt = 0.001;
	double T = ntime;
	double tau = 8;
	double Di = 1;
	int N = int(T/dt);
	double x[N];
	double x2[N];
	

	double tau2 = 0.1;
	

	vector <double> * Neurons= new vector <double> [nunits];


	ofstream fout("output.txt");
	double *diff;
	diff = OrnsteinUhlebeck(x2,dt, T, tau, Di);
	fout << 0;
	for (int i = 0; i < nunits; i ++){
		
		double * poi;
		poi = OrnsteinUhlebeck2(x,dt, T, tau2, diff);


		for (int s = 0; s < N; s ++){


			Neurons[i].push_back(*(poi+s));
			//cout <<(*(poi + s))<< endl;
			//SaveFile << " ";
			fout << ", "<< *(poi+s);
			//fout << " ";
			
		}
		
//fout << "\n";


}




fout.close();
}



double* OrnsteinUhlebeck2(double *x,double dt, double T, double tau, double* D) {

	int N = int(T/dt);
	x[0] = 0;
	for (int t = 0; t < N; t ++){
		x[t+1] = x[t] -x[t]*dt/tau + TMath::Sqrt(*(D+t)*dt)*gRandom->Gaus(0,1);
	}

	return x;
	
	
}

double* OrnsteinUhlebeck(double *x, double dt, double T, double tau, double D) {

	int N = int(T/dt);
	x[0] = 0;
	for (int t = 0; t < N; t ++){
		x[t+1] = x[t] -x[t]*dt/tau + TMath::Sqrt(D*dt)*gRandom->Gaus(0,1);
	}


	for (int i = 0; i < N; i ++){
		if (x[i] < 0.02){
			x[i] = 0.02;
		}
	}
	return x;
	
	
}


    
	


	




