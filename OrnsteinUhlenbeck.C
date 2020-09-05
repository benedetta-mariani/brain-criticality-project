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


vector <double> OUGillespieDiffusion(double dt, double T, double tau, double D);


void Main(int nunits, int T){ 
	double dt = 0.001;
	double tauD = 6;
	double Di = 1;
	double tau = 0.08;
	int N = int(T/dt);
	int n = int(tauD/dt*100) + N; 
	
	vector <double> diff;
	diff = OUGillespieDiffusion(dt, T, tauD, Di);
	
	ofstream fout("output.csv");

	//////////////////////////////////////////ADJACENCY MATRIX//////////////////////////////////////////
	/// If you want to change the adjacency matrix make sure that its spectral radius is less than 1 ///

	double a[nunits][nunits];

	for (int r = 0; r < nunits; r++){
	    for (int s = 0; s < r; s++){
	        double g = gRandom->Uniform(0,1);
	        if (g <= 0.5){
	            a[r][s] = 1;
	        }
	        else{
	            a[r][s] = 0;
	        }   
	        a[s][r] = a[r][s];
	     }    
	    a[r][r] = 0; // no self edges
	 }


	for (int r = 0; r < nunits; r++){
	    for (int s = 0; s < r + 1; s++){
	        a[r][s] = a[r][s]/120; //Divided by this constant so to have spectral radius of a < 1. 
	        a[s][r] = a[r][s];
	    }
	}    
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	vector <double> x[nunits];

	for (int i = 0; i < nunits; i ++){

		x[i].push_back(0.);
	}
		
	double number;
	for (int t = 0; t < N-1; t ++){

		for (int i = 0; i < nunits; i ++){
			double value = 0;
			for (int j = 0; j < nunits; j ++){
				value += a[i][j]*x[j][t]*dt/tau;
			}

			number = x[i][t] - x[i][t]*dt/tau + value + TMath::Sqrt(diff[int(tauD/dt*100)+t]*dt)*gRandom->Gaus(0,1);
			x[i].push_back(number);
		}

	}



	for (int i = 0; i < nunits; i ++){
		fout << 0;
		for (int j = 1; j < N; j ++){
			fout << "," << x[i][j];
		}
		fout << "\n";
	}
		


	fout.close();
	
}



vector <double> OUGillespieDiffusion(double dt, double T, double tau, double D){
	
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
		if (x[i] < 0.2){
			x[i] = 0.2;
		}
	}	
	return x;

}

		



    
	


	




