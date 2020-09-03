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


double* OUGillespieDiffusion(double dt, double T, double tau, double D);
double* OrnsteinUhlebeck(double dt, double T, double tau, double* D);


void Main(int nunits, int T){ 

	double dt = 0.001;
	double tau = 10;
	double Di = 1;
	double tau2 = 0.08;
	int N = int(T/dt);

	double* diff;
	diff = OUGillespieDiffusion(dt, T, tau, Di);


	ofstream fout("output.csv");
	
	for (int i = 0; i < nunits; i ++){
		
		double* poi;
		poi = OrnsteinUhlebeck(dt, T, tau2, diff);
		fout << *(poi+ int(100*tau2/dt));  // remove initial transient

		for (int s = int(100*tau2/dt)+1; s < N; s ++){
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
	for (int t = 0; t < N-1; t ++){
		x[t+1] = x[t] -x[t]*dt/tau + TMath::Sqrt(*(D+int(tau/dt*100)+t)*dt)*gRandom->Gaus(0,1);
	}

	return x;
	
	
}

double* OUGillespieDiffusion(double dt, double T, double tau, double D){


	int N = int(T/dt);
	int n = int(tau/dt*100) + N; // to avoid initial transient
	double *x = new double[n];
	x[0] = 0;
    double mu = TMath::Exp(-dt/tau);
    double sigma = TMath::Sqrt(D*tau/2*(1-pow(mu,2)));

    for (int t = 0; t < n-1; t ++){
        x[t+1] = x[t]*mu + sigma*gRandom->Gaus(0,1);
     }
  
	for (int i = 0; i < n; i ++){
		if (x[i] < 0.2){
			x[i] = 0.2;
		}
	}
	
return x;

}

	
	



    
	


	




