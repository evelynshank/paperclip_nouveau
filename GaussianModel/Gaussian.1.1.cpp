// SphHarEvol.1.cpp.cpp : Defines the entry point for the console application.
// Written by: Suren Gourapura
// Date: 10/28/17
// Goal: The goal is to take a given direction and standard deviation and to create a gain pattern

#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
using namespace std;

double distSph(int theta, int phi, double thetaM, double phiM) {
	const double PI = 3.14159;
	double distSqX = pow(sin(theta*PI / 180.0)*cos(phi*PI / 180.0) - sin(thetaM*PI / 180.0)*cos(phiM*PI / 180.0), 2);
	double distSqY = pow(sin(theta*PI / 180.0)*sin(phi*PI / 180.0) - sin(thetaM*PI / 180.0)*sin(phiM*PI / 180.0), 2);
	double distSqZ = pow(cos(theta*PI / 180.0) - cos(thetaM*PI / 180.0), 2);
	return sqrt(distSqX + distSqY + distSqZ);
}

double normPdf(double x, double mu, double sigma){
	const double PI = 3.14159;
	return exp(-1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * PI));
}

int main()
{
double thetaM; // M for mean value
double phiM;
double stdev;
const double PI = 3.14159;
double norm = 0;
double gain[36][72];

// Enter the direction and standard deviation of the signal
cout << "Enter the direction and standard deviation (all in degrees) in the following format: theta phi stdev" << endl;
cin >> thetaM;
cin >> phiM,
cin >> stdev;

// We take each point on the sphere, find out its distance from (thetaM, phiM), and assign it a gain value based on the value of a normal PDF at that distance

for (int i = 0; i < 36; i++) {
	for (int j = 0; j < 72; j++) {
		double dist = distSph(i*5, j*5, thetaM, phiM) / 2;
		gain[i][j] = normPdf(dist, 0, stdev/ 72);
		//gaindB[i][j] = 10 * log10( normPdf(dist, 0, stdev/360) / ( norm /(360 * 180) ));
	}
}

// Now we normalize it in two steps, first calculate the sum over all points and then divide all value so that the sum  over all points equals 2592 (2592 = 36 * 72)

for (int i = 0; i < 36; i++) {
	for (int j = 0; j < 72; j++) {
		norm += gain[i][j] * sin(5* i * PI / 180.0)* 2 * PI; // We mulitply by the circumference at each point to normalize correctly
	}
}
cout << norm << endl;

for (int i = 0; i < 36; i++) {
	for (int j = 0; j < 72; j++) {
		gain[i][j] = 2592 * gain[i][j] / norm;
	}
}

// Finally, we output these values in a Mathematica-Friendly format 

for (int i = 0; i < 36; i++) {
	for (int j = 0; j < 72; j++) {
		cout << setprecision(10) << fixed <<"{" << i * 5 * PI / 180.0  << ", " << j * 5 * PI / 180.0 << ", " << gain[i][j] << "}, ";
	}
}

return 0;
}
