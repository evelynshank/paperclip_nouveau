#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main(){
	int numElements = 2664;
	int numCoeff = 13;
	ofstream myfile;
	myfile.open("FitTest.P83.33.Y00-Y012.txt");
	double theta[numElements];
	double phi[numElements];
	double r[numElements];
	double sphHar[numElements];
	double Pi = M_PI;
	double AvgDev=0.0;
	double Chi2=0.0;
	double coeff[numCoeff];
	
// Obtain the Coefficients
	cout << "Enter the coefficients for L=0 L=1 .... L=12 with enters in between:"<< endl;
	fill(coeff, coeff+numCoeff, 0.0);
	for(int i=0; i< numCoeff; i++){
		cin >> coeff[i];
	}
// Obtain the r (Gain or Phase) data
	cout << "Enter r: "<< endl;
	fill(r, r+numElements, 0.0);
	for(int i=0; i< numElements; i++){
		cin >> r[i];
	}

// Populate the theta array
	fill(theta, theta+numElements, 0.0);
	for (int j = 0; j < 72; j++){
		for (int k=0; k<37; k++){
			theta[j*37+k]=5*k*Pi/180.0;
		}
	}

// Populate the phi array
	fill(phi, phi+numElements, 0.0);
	for (int l = 0; l < 72; l++){
		for (int m=0; m < 37; m++){
			phi[l*37+m]=5*l*Pi/180.0;
		}
	}
	
// Populate harmonic data points
	fill(sphHar, sphHar+numElements, 0.0);
	for (int p = 0; p < numElements; p++){
		sphHar[p]= 
//Y00		
		coeff[0]*(1/2.0)*(1/sqrt(Pi)) 
//Y01
		+coeff[1]*pow(10,-6)*(1/2.0)*sqrt(3/Pi)*cos(theta[p])
//Y02
		+coeff[2]*(1/4.0)*sqrt(5/Pi)*(3*pow(cos(theta[p]), 2)- 1)
//Y03
		+coeff[3]*pow(10, -7)*(1/4.0)*sqrt(7/Pi)*(5*pow(cos(theta[p]),3)- 3*cos(theta[p]))
//Y04
		+coeff[4]*(3/16.0)*sqrt(1/Pi)*(35*pow(cos(theta[p]),4) - 30*pow(cos(theta[p]),2)+3)
//Y05
		+coeff[5]*pow(10,-7)*(1/16.0)*sqrt(11/Pi)*(15*cos(theta[p]) - 70*pow(cos(theta[p]),3)+63*pow(cos(theta[p]),5))
//Y06
		+coeff[6]*(1/32.0)*sqrt(13/Pi)*(-5 + 105*pow(cos(theta[p]),2)-315*pow(cos(theta[p]),4) + 231*pow(cos(theta[p]),6))
//Y07					
		+coeff[7]*pow(10, -7)*(1/32.0)*sqrt(15/Pi)*(-35*cos(theta[p])+ 315*pow(cos(theta[p]),3) -693*pow(cos(theta[p]),5) + 429*pow(cos(theta[p]),7))
//Y08					
		+coeff[8]*(1/256.0)*sqrt(17/Pi)*(35 - 1260*pow(cos(theta[p]),2) + 6930*pow(cos(theta[p]),4) - 12012*pow(cos(theta[p]),6) + 6435*pow((cos(theta[p])),8))
//Y09
		+coeff[9]*(1/256.0)*sqrt(19/Pi)*(315*cos(theta[p])- 4620*pow(cos(theta[p]),3) + 18018*pow(cos(theta[p]),5) - 25740*pow(cos(theta[p]),7) + 12155*pow((cos(theta[p])),9))
//Y10					
		+coeff[10]*(1/512.0)*sqrt(21/Pi)*(-63 +3465*pow(cos(theta[p]),2) - 30030*pow(cos(theta[p]),4) + 90090*pow(cos(theta[p]),6) -109395*pow((cos(theta[p])),8)+46189*pow(cos(theta[p]),10))		
//Y011
		+coeff[11]*(1/512.0)*sqrt(23/Pi)*(-693*pow(cos(theta[p]),1) +15015*pow(cos(theta[p]),3) - 90090*pow(cos(theta[p]),5) +218790*pow((cos(theta[p])),7)-230945*pow(cos(theta[p]),9)+88179*pow(cos(theta[p]),11))
//Y012
		+coeff[12]*(1/2048.0)*sqrt(25/Pi)*(231 -18018*pow(cos(theta[p]),2) +225225*pow(cos(theta[p]),4) - 1021020*pow(cos(theta[p]),6) +2078505*pow((cos(theta[p])),8)-1939938*pow(cos(theta[p]),10)+676039*pow(cos(theta[p]),12))
		;
	}

// Calculate AvgDev
	for(int n=0; n < numElements; n++){
		AvgDev += abs(r[n]-sphHar[n]);
	}
	AvgDev = AvgDev/numElements;
	
// Calculate Chi^2
	for(int n=0; n < numElements; n++){
		Chi2 += pow((r[n]-sphHar[n]),2)/sphHar[n];
	}
// Close file
	myfile << "Average Deviation (Delta R / N): "<< AvgDev << endl;
	myfile << "Chi^2: "<< Chi2 << endl;
	myfile.close();
	return 0;
}
