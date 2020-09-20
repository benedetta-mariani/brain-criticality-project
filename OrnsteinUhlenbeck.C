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


vector <double> OUGillespieDiffusion(double dt, double T, double tau, double D, double Dstar);
vector<vector<double>> AdjacencyMatrix(double pstar, int nunits);
int MinimalCheck(vector<vector<double>> a, int nunits);
vector<vector<double>> ConfigModel(int nunits);
vector<vector<double>> units(int nunits, vector<vector<double>> a, vector<double> diff,double tauD, double tau, double dt, bool connected, int N);


void Main(int nunits, int T){ 
	double dt = 0.001;
	double tauD = 6;
	double Di = 1;
	double tau = 0.08;
	double Dstar = 0.2;
	double pstar = 0.5;
	bool connected = kFALSE;
	int N = int(T/dt);
	int n = int(tauD/dt*100) + N; 
	
	vector <double> diff;
	diff = OUGillespieDiffusion(dt, T, tauD, Di, Dstar);
	vector<vector<double>> a(nunits, vector<double> (nunits,0.));
	//a = AdjacencyMatrix(pstar, nunits);
	a = ConfigModel(nunits);
	if (MinimalCheck(a, nunits) > 0){
		int f = 1;
		do {
			for (int i = 0; i < nunits; i ++){
				for (int j = 0; j < i+1; j++){
					a[i][j] = a[i][j]/(2*f);
					a[j][i] = a[i][j];
					f ++; 
				}
			}
		}
		while(MinimalCheck(a, nunits) > 0);
	}

	ofstream fout("output.csv");
	vector<vector<double>> neuro(nunits, vector<double> (N,0.));
	neuro = units(nunits, a, diff,tauD,tau, dt, connected, N);
	
	//diff[int(tauD/dt*100)
	for (int i = 0; i < nunits; i ++){
		fout << neuro[i][0];
		for (int j = 1; j < N; j ++){
			fout << "," << neuro[i][j];
		}
		fout << "\n";
	}
	fout.close();
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


vector<vector<double>> AdjacencyMatrix(double pstar, int nunits){


	vector<vector<double>> a(nunits, vector<double> (nunits, 0.));

	for (int r = 0; r < nunits; r++){
	    for (int s = 0; s < r; s++){
	        double g = gRandom->Uniform(0,1);
	        if (g <= pstar){
	            a[r][s] = 1;
	        }
	        else{
	            a[r][s] = 0;
	        }   
	        a[s][r] = a[r][s];
	     }    
	    a[r][r] = 0;
	 }


	
	for (int r = 0; r < nunits; r++){
	    for (int s = 0; s < r + 1; s++){
	        a[r][s] = a[r][s]/120; //Divided by this constant so to have spectral radius of a < 1. 
	        a[s][r] = a[r][s];
	    }
	}    



	return a;
}


int MinimalCheck(vector<vector<double>> a, int nunits){
	double m = 0; 
	double k  =  static_cast<double> (nunits);
	
	for (int r = 0; r < nunits; r++){
	    for (int s = 0; s < nunits; s++){
	    	m += a[r][s];
	    }
	} 
	m = m/2;
	double bound = TMath::Sqrt( 2*m - nunits - k + 5/2 + TMath::Sqrt(2*m - 2* nunits + 9/4));
	if (bound> 1){
		cout << "There may be problems" << endl;
		cout << "The upper bound value is..." << bound << endl;
		return 1;
	}
	else{
		cout << "Surely there aren't problems" << endl;
		cout << "The upper bound value is..." << bound << endl;
		return 0;
	}
}

		

vector<vector<double>> ConfigModel(int nunits){

	vector<int> node_degree(nunits, 0);


	int i,j;
    int n_links;
    int max_size; //Maximum size if we want an uncorrelated network
    double gamma = 2.5;

    vector<int> link_vector;

    max_size = sqrt(nunits); //Max size to avoid correlatons

    n_links = 0; //There's no link yet

    //Compute all the quantities we need to generate the degrees,
    double kmax = pow(max_size, 1.0-gamma);
    double kmin = pow(1, 1.0-gamma);
    double invgamma = 1.0 / (1.0 - gamma);

    double ran;
    for (i=0; i < nunits; i++){
    	ran = gRandom-> Uniform(0,1);
        node_degree[i] = floor( pow( ran*(kmax - kmin) + kmin, invgamma ) ); //Generate degrees
        n_links += node_degree[i]; //Update number of links
    }
  
    //Make sure we have even number of links
    if (n_links % 2 == 1)
    {
        node_degree[0] += 1;
        n_links += 1;
    }

	
	
	vector<vector<double>> a(nunits, vector<double> (nunits, 0.));

	for (int i = 0; i < nunits; i ++){
		for (int j = 0; j < i +1; j ++){
			ran = gRandom->Uniform(0,1);
			if (ran <= (node_degree[i]*node_degree[j])/(2*n_links -1) ){
				a[i][j] = 1;
			}
			else{
				a[i][j] = 0;
			}
			a[j][i] = a[i][j];
		}
	}

	return a;

}

vector<vector<double>> units(int nunits, vector<vector<double>> a, vector<double> diff, double tauD, double tau, double dt, bool connected, int N){
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
					value += a[i][j]*x[j][t]*dt/tau;
				}
			}
			x[i][t+1] = x[i][t] - x[i][t]*dt/tau + value + TMath::Sqrt(diff[tauD*100/dt+ t]*dt)*gRandom->Gaus(0,1);
			
		}
	}
	//cout << x[0].size() << endl;
	return x;

}



    
	


	




