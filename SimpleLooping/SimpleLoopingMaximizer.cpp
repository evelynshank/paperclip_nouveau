
// Written by: Suren Gourapura
// Date: 3/12/18
/* Goal: To see how well the SphHarEvolRing.1.3.cpp code works, we write this competitor code whose
		job is to achieve a given fitness score for a given theta and phi by repeadedly randomly 
		generating 13 coefficients until one combination gets a score higher than asked for.
*/
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
#include <vector>
#include <fstream>
using namespace std;

double Pi = 3.14159265359;
int SphHarMAX = 13;

	
double integrateSpecies( double species[], double Theta, double Spread){
	
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
			species[0]*(1/2.0)*(1/sqrt(Pi)) 
//Y01
			+species[1]*(1/2.0)*sqrt(3/Pi)*cos(degree)
//Y02
			+species[2]*(1/4.0)*sqrt(5/Pi)*(3*pow(cos(degree), 2) - 1)
//Y03
			+species[3]*(1/4.0)*sqrt(7/Pi)*(5*pow(cos(degree),3) - 3*cos(degree))
//Y04
			+species[4]*(3/16.0)*sqrt(1/Pi)*(35*pow(cos(degree),4) - 30*pow(cos(degree),2) + 3)
//Y05
			+species[5]*(1/16.0)*sqrt(11/Pi)*(63*pow(cos(degree),5) - 70*pow(cos(degree),3) + 15*cos(degree))
//Y06
			+species[6]*(1/32.0)*sqrt(13/Pi)*(231*pow(cos(degree),6) - 315*pow(cos(degree),4) + 105*pow(cos(degree),2) - 5 )
//Y07					
			+species[7]*(1/32.0)*sqrt(15/Pi)*(429*pow(cos(degree),7) - 693*pow(cos(degree),5) + 315*pow(cos(degree),3) - 35*cos(degree))
//Y08					
			+species[8]*(1/256.0)*sqrt(17/Pi)*(6435*pow((cos(degree)),8) - 12012*pow(cos(degree),6) + 6930*pow(cos(degree),4) - 1260*pow(cos(degree),2) + 35)
//Y09
			+species[9]*(1/256.0)*sqrt(19/Pi)*(12155*pow((cos(degree)),9) - 25740*pow(cos(degree),7) + 18018*pow(cos(degree),5) - 4620*pow(cos(degree),3) + 315*cos(degree))
//Y10					
			+species[10]*(1/512.0)*sqrt(21/Pi)*(46189*pow(cos(degree),10) - 109395*pow((cos(degree)),8) + 90090*pow(cos(degree),6) - 30030*pow(cos(degree),4) + 3465*pow(cos(degree),2) - 63)		
//Y011
			+species[11]*(1/512.0)*sqrt(23/Pi)*(88179*pow(cos(degree),11) - 230945*pow(cos(degree),9) + 218790*pow((cos(degree)),7) - 90090*pow(cos(degree),5) + 15015*pow(cos(degree),3) - 693*cos(degree))
//Y012
			+species[12]*(5/2048.0)*sqrt(1/Pi)*(676039*pow(cos(degree),12) - 1939938*pow(cos(degree),10) + 2078505*pow((cos(degree)),8) - 1021020*pow(cos(degree),6) + 225225*pow(cos(degree),4) - 18018*pow(cos(degree),2) + 231)
			, 2)*sin(degree);
	}
// Finally, we give the summed values of the spherical harmonic from theta min to theta max as the test score (a rough kind of integration by .5 degree increments)			
	return Sum;	
}


int main()
{
	double species[SphHarMAX], bestSpecies[SphHarMAX];
	double Theta, Spread;
	double goalVal;
	double bestVal=0;
	int iteration = 0;
	srand(time(NULL));
	ofstream myfile;
	myfile.open("SimpleLooping_0.2_3.12.18.txt");
/*
// Enter the goal theta and spread
	cout << "Enter the goal theta and spread (in degrees), and goal fitness value in the following format: theta spread fitness" << endl;
	cin >> Theta;
	cin >> Spread;
	cin >> goalVal;
*/
	Theta = 90;
	Spread = 20;
	goalVal = 0.775;	
	while(bestVal < goalVal){
		iteration++;
		for (int j = 0; j <= SphHarMAX - 1; j++){
			species[j] = ((double)rand() / (double)(RAND_MAX))*2 - 1; // gives a random value between -1 and 1
		}
		
		double counter = pow(integrateSpecies(species, 90, 180)*0.0548594, 0.5); // Normalize the randomly generated species
		for (int j = 0; j <= SphHarMAX - 1; j++){
			species[j] = species[j] / counter;
		}
		
		if (integrateSpecies(species, Theta, Spread)*0.0548594 > bestVal)	{
			
			bestVal = integrateSpecies(species, Theta, Spread)*0.0548594; // Assign a new best score
		
			for (int j = 0; j <= SphHarMAX - 1; j++){
				bestSpecies[j] = species[j];				// Assign a new best species
			}
		}
		if (iteration%1000 == 0){
			cout << "Iteration #: " << iteration << endl;
		}
	}
		
// The work is done, now time to see the results of the final generation! We follow the same ranking protocol
	myfile << "The best Fitness Value was: " << bestVal << endl;
	myfile << "The final iteration was: " << iteration << endl;
	myfile << "The Species coefficients are: "<< endl;

	for (int j = 0; j <= SphHarMAX-1; j++){
		myfile << "m1Y0" << j << " = " << bestSpecies[j] << endl;
	}
	myfile.close();

return 0;
}
