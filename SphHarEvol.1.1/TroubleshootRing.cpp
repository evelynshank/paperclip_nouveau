
// Written by: Suren Gourapura
// Date: 10/28/17
// Goal: The goal is to take a given theta angle and spread in theta to evolve an antenna using spherical harmonics that best satisfies it.

#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
#include <vector>
using namespace std;

double Pi = 3.14159265359;
int SphHarMAX = 13;
int PopMAX = 100;


// Code modified from https://codereview.stackexchange.com/questions/110793/insertion-sort-in-c
void insertionSort(double array[], int length){
    int i,j;
    for (i = 1; i < length; i++) {
        double temp = array[i];
        for (j = i; j > 0 && array[j - 1] < temp; j--) {
            array[j] = array[j - 1];
        }
        array[j] = temp;
    }
}


double integrateSpecies( vector<vector<double> >& pop, int species, double Theta, double Spread){
	
	double Min, Max;
// Make sure that the maximum is less than or equal to 180 degrees and the minimum is greater than or equal to 0
	if (Theta + Spread/2.0 > 180){Max = 180;}
	else{Max = Theta + Spread/2.0;}
			
	if (Theta - Spread/2.0 < 0){Min = 0;}
	else{Min = Theta - Spread/2.0;}
	
	double Sum = 0;			
// Sum the spherical harmonic values by half-degree increments from the Min to Max theta values
	for (double degree = Min * Pi / 180; degree <= Max *Pi / 180; degree = degree + Pi/360){
		Sum += pow(
//Y00		
			pop[species][0]*(1/2.0)*(1/sqrt(Pi)) 
//Y01
			+pop[species][1]*(1/2.0)*sqrt(3/Pi)*cos(degree)
//Y02
			+pop[species][2]*(1/4.0)*sqrt(5/Pi)*(3*pow(cos(degree), 2) - 1)
//Y03
			+pop[species][3]*(1/4.0)*sqrt(7/Pi)*(5*pow(cos(degree),3) - 3*cos(degree))
//Y04
			+pop[species][4]*(3/16.0)*sqrt(1/Pi)*(35*pow(cos(degree),4) - 30*pow(cos(degree),2) + 3)
//Y05
			+pop[species][5]*(1/16.0)*sqrt(11/Pi)*(63*pow(cos(degree),5) - 70*pow(cos(degree),3) + 15*cos(degree))
//Y06
			+pop[species][6]*(1/32.0)*sqrt(13/Pi)*(231*pow(cos(degree),6) - 315*pow(cos(degree),4) + 105*pow(cos(degree),2) - 5 )
//Y07					
			+pop[species][7]*(1/32.0)*sqrt(15/Pi)*(429*pow(cos(degree),7) - 693*pow(cos(degree),5) + 315*pow(cos(degree),3) - 35*cos(degree))
//Y08					
			+pop[species][8]*(1/256.0)*sqrt(17/Pi)*(6435*pow((cos(degree)),8) - 12012*pow(cos(degree),6) + 6930*pow(cos(degree),4) - 1260*pow(cos(degree),2) + 35)
//Y09
			+pop[species][9]*(1/256.0)*sqrt(19/Pi)*(12155*pow((cos(degree)),9) - 25740*pow(cos(degree),7) + 18018*pow(cos(degree),5) - 4620*pow(cos(degree),3) + 315*cos(degree))
//Y10					
			+pop[species][10]*(1/512.0)*sqrt(21/Pi)*(46189*pow(cos(degree),10) - 109395*pow((cos(degree)),8) + 90090*pow(cos(degree),6) - 30030*pow(cos(degree),4) + 3465*pow(cos(degree),2) - 63)		
//Y011
			+pop[species][11]*(1/512.0)*sqrt(23/Pi)*(88179*pow(cos(degree),11) - 230945*pow(cos(degree),9) + 218790*pow((cos(degree)),7) - 90090*pow(cos(degree),5) + 15015*pow(cos(degree),3) - 693*cos(degree))
//Y012
			+pop[species][12]*(5/2048.0)*sqrt(1/Pi)*(676039*pow(cos(degree),12) - 1939938*pow(cos(degree),10) + 2078505*pow((cos(degree)),8) - 1021020*pow(cos(degree),6) + 225225*pow(cos(degree),4) - 18018*pow(cos(degree),2) + 231)
			, 2)*sin(degree);
	}
// Finally, we give the summed values of the spherical harmonic from theta min to theta max as the test score (a rough kind of integration by .5 degree increments)			
	return Sum;	
}


int main()
{
	double pop[PopMAX][SphHarMAX];
	double Theta, Spread;
	int Gen = 0;
	srand(time(NULL));
	
// Enter the goal theta and spread
	cout << "Enter the goal theta and spread (in degrees), and generations in the following format: theta spread generations" << endl;
	cin >> Theta;
	cin >> Spread;
	cin >> Gen;
	
//  fill the population
	for (int i = 0; i <= PopMAX - 1; i++){
	//	pop[i][0] = 1;
	//	pop[i][1] = 2;
		for (int j = 0; j <= SphHarMAX - 1; j++){
			pop[i][j] = i*j; 
		}
	}

	for (int i = 0; i <= PopMAX - 1; i++){
		cout << i << " : ";
		for (int j = 0; j <= SphHarMAX - 1; j++){
			cout << pop[i][j] << " "; 
		}
		cout << endl;
	}
	
	for (int i = 1; i <= 1; i++){
		for (int j = 0; j <= SphHarMAX-1; j++){
			cout << "m1Y0" << j << " = " << pop[i][j] << endl;
		}
		cout << endl;
	}	
	

	vector<vector<double> >vPop;
	vPop.resize(PopMAX);
	for(int i=0; i<PopMAX; i++){
		vPop[i].resize(SphHarMAX);
	}
	for(int i=0; i<PopMAX; i++){
		for(int j=0; j<SphHarMAX; j++){
			vPop[i][j] = pop[i][j];
		}
	}
// Now we normalize	
	for (int i = 0; i <= PopMAX - 1; i++){
		double counter = pow(integrateSpecies(vPop, i, 90, 180)*0.0548594, 0.5);
		for (int j = 0; j <= SphHarMAX - 1; j++){
			pop[i][j] = pop[i][j] / counter;
		}		
		cout << i << " counter: " << counter << endl;
	}

	for (int i = 0; i <= PopMAX - 1; i++){
		cout << i << " : ";
		for (int j = 0; j <= SphHarMAX - 1; j++){
			cout << pop[i][j] << " "; 
		}
		cout << endl;
	}

	for (int i = 1; i <= 2; i++){
		for (int j = 0; j <= SphHarMAX-1; j++){
			cout << "m1Y0" << j << " = " << pop[i][j] << endl;
		}
		cout << endl;
	}


return 0;
}
/*
// Evolution Algorithim 4: Introduce 10 random species into the population
// First step is to give each value in each species random number between -1 and 1. Then, we normalize the species values to 1
		for (int i = 0; i <= 9; i++){
			for (int j = 0; j <= SphHarMAX - 1; j++){
				nextPop[90 + i][j] = ((double)rand() / (double)(RAND_MAX))*2- 1; // gives a random value between -1 and 1
			}
		}		
		
// Now we normalize nextPop: We divide all coefficients by their integral over all space to normalize to 1
// To pass this population off to functions, we make a 2D vector with the same values called vPop
		for(int i=0; i<PopMAX; i++){
			for(int j=0; j<SphHarMAX; j++){
				vPop[i][j] = nextPop[i][j];
			}
		}
// Now we normalize	
		for (int i = 0; i <= PopMAX - 1; i++){
			double counter = integrateSpecies(vPop, i, 90, 180);
			for (int j = 0; j <= SphHarMAX - 1; j++){
				nextPop[i][j] = nextPop[i][j] / counter;
			}		
		}
		*/
