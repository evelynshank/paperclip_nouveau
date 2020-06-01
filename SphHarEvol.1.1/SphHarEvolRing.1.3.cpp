
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
	
/*// Enter the goal theta and spread
	cout << "Enter the goal theta and spread (in degrees), and generations in the following format: theta spread generations" << endl;
	cin >> Theta;
	cin >> Spread;
	cin >> Gen;
*/
	Theta = 90;
	Spread = 20;
	Gen = 100;
// Randomly fill the population. First step is to give each value in each species random number between -1 and 1. Then, we normalize the species values to 1
	for (int i = 0; i <= PopMAX - 1; i++){
		for (int j = 0; j <= SphHarMAX - 1; j++){
			pop[i][j] = ((double)rand() / (double)(RAND_MAX))*2 - 1; // gives a random value between -1 and 1
		}
	}
	
// We divide all coefficients by their integral over all space to normalize to 1
// To pass this population off to functions, we make a 2D vector with the same values called vPop
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
	}



// Now, we evolve. We loop the following code for each generation
	for (int g = 1; g <= Gen; g++){
	//	cout << endl << "Generation: " << g << endl;
		
// The first objective is to calculate how well our current population's species are doing
		double testScores[100];
// To pass this population off to functions, we again make a 2D vector with the same values called vPop
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

// Now, to test the curent population, we integrate each species from [theta-spread/2] to [theta + spread/2]
		for (int i = 0; i <= PopMAX - 1; i++){
			testScores[i] = integrateSpecies(vPop, i, Theta, Spread)*0.0548594;
		}

// The following matrix will hold the next generation's species. In the end of the evolution, we will make pop = nextPop, so the next evolution afterwards will act on nextPop 
		double nextPop[PopMAX][SphHarMAX] = {};

		
// Evolution Algorithim 1: Take the 10 species with the best scores and pass them onto nextPop

// We want to organize the population by their testScores, so we make a temporary matrix to hold the ranked initial population:rankedPop 
// Also, a temporary ordered testScores array: rankedTestScores
		double rankedPop[PopMAX][SphHarMAX] = {};
		double rankedTestScores[PopMAX] = {};
// First, equate the test scores with the ranked test scores. Then we sort the ranked test scores: greatest to least
		for (int i=0; i <= PopMAX-1; i++){
			rankedTestScores[i] = testScores[i];
		}

		insertionSort(rankedTestScores, PopMAX);
		
// Next, we find where testScores = rankedTestScores to place the pop species in the right place in rankedPop
		for (int i=0; i <=PopMAX-1; i++){
			for (int k=0; k <=PopMAX-1; k++){
				if (testScores[k] == rankedTestScores[i]){
					
					for (int l=0; l <=SphHarMAX-1; l++){ // We need to copy the whole of the species in the kth position to rankedPop in the ith position
						rankedPop[i][l] = pop[k][l];
					}
					
				}
			}
		}
	
// We print out the highest ranking species's score and it's array in mathematica format, for easy plotting	
//		cout <<rankedTestScores[0]<< endl;	
//		cout << "Generation: " << g << " First Score in Mathematica Format with score "<< rankedTestScores[0]<< " : " << endl; 
//		for (int i = 0; i <= SphHarMAX-1; i++){
//			cout << "m1Y0" << i<< " = " << rankedPop[0][i] << endl;
//		}*/
// Finally, we copy over the top 10 best species from pop to nextPop		
		for (int i = 0; i <= 10-1; i++){
			for (int j = 0; j <= SphHarMAX-1; j++){
				nextPop[i][j]= rankedPop[i][j];
			}

		}



// Evolution Algorithim 2: Take 10 random species, find the one with the best score, and [randomly mutate one of it's array values to obtain an offspring] 10 times. 
// Do this whole algorithim 3 times
		for(int a2 = 1; a2 <= 3; a2++){		
		
			int choose10[10]; // Create an array with 10 random values, 0-99. This array determines the 10 random species that will undergo a tournament selection (a.k.a. simply choosing the highest score species out of the 10)
			for (int i = 0; i <= 9; i++){
				choose10[i]= rand() % 100;
			}
		
// Since we have pop already organized from best to worst, we simply find the lowest value in choose10 to find the winner of the tournament. The winning species is called: Algorithim 2 Best Value
			int Alg2BestVal = choose10[0]; // Assume the best species is the first one
		
			for (int i = 1; i <= 9; i++){
				if (Alg2BestVal > choose10[i]){  // If another species in the choose10 array has a lower value, it is now the best species
					Alg2BestVal = choose10[i];		
				}
			}
			
// Now, we mutate one part of the species's array, pop[Alg2BestVal][], to create an offspring. We do this 10 times, normalizing after each one
			int mutateLocation[10]; // Create the locations for mutation of the 10 offspring
			for (int i = 0; i <= 9; i++){
				mutateLocation[i]= rand() % SphHarMAX;
			}
			
// We do the process below 10 times, (to spots 10*a2 to 19*a2 in nextPop)
			for (int i = 0; i <= 9; i++){
			
				for (int j = 0; j <= SphHarMAX-1; j++){
					nextPop[i + 10*a2][j]= rankedPop[Alg2BestVal][j];	// Copy the Alg2BestVal species over. Note: a2 is added to be able to do the whole of Algorithm 2, 3 times
				}
				for (int j = 0; j <= SphHarMAX-1; j++){
					if(mutateLocation[i] == j){  // make sure that the location to be mutated is chosen randomly, by using mutateLocation[]
						nextPop[i + 10*a2][j] = ((double)rand() / (double)(RAND_MAX))-.5;// Mutate this location	
					}
				}	
			}
// We will normalize the whole next population after all evolutionary algorithims			
		}
		
		
// Evolution Algorithim 3: Take 20 random species, run two seperate tournaments to find 2 parents, and [swap a random array location with each other to obtain two offspring (for both combinations)] 5 times.
// We do this Algorithm 5 times, for a total of 50 offspring		
		for(int a3 = 1; a3 <= 5; a3++){		
			
			int choose10A[10], choose10B[10]; // Create 2 arrays with 10 random values, 0-99. These arrays determine the 20 random species that will undergo two seperate tournament selections
			int Alg3BestValA = 0, Alg3BestValB = 0;
			
			while (Alg3BestValA == Alg3BestValB){ // To make sure that the 2 parents aren't the same species, we run the following as long as they are the same
				for (int i = 0; i <= 9; i++){
					choose10A[i]= rand() % 100;
					choose10B[i]= rand() % 100;
				}
		
// Since we have pop already organized from best to worst, we simply find the lowest value in choose10 to find the winners of the tournaments. The winning species are called: Algorithim 2 Best Value A/B
				Alg3BestValA = choose10A[0]; // Assume the best species are the first ones
				Alg3BestValB = choose10B[0];
				for (int i = 1; i <= 9; i++){
					if (Alg3BestValA > choose10A[i]){  // If a lower value is found, make best value that lower value
						Alg3BestValA = choose10A[i];		
					}
					if (Alg3BestValB > choose10B[i]){ 
						Alg3BestValB = choose10B[i];		
					}
				}					
			}
// Now, we swap one part of the species's arrays in both ways to create 2 offspring. We do this 5 times, normalizing after each one
			int swapLocation[5]; // Create the locations for swapping
			for (int i = 0; i <= 4; i++){
				swapLocation[i]= rand() % SphHarMAX;
			}
	
// We do the process below 5 times, (to spots 40*a3 to 49*a3 in nextPop)
			for (int i = 0; i <= 4; i++){
// First, we copy over the two parents in the 10 offspring spots in nextPop, in an A,B,A,B,... pattern			
				for (int j = 0; j <= SphHarMAX-1; j++){
					nextPop[(i*2) + 10*a3 + 30][j]= rankedPop[Alg3BestValA][j];	// Copy the Alg3BestVal[A and B] species over. Note: a3 is added to be able to do the whole of Algorithm 3, 5 times
					nextPop[(i*2 + 1) + 10*a3 + 30][j]= rankedPop[Alg3BestValB][j];
				}
				for (int j = 0; j <= SphHarMAX-1; j++){
					if(swapLocation[i] == j){  // make sure that the location to be swapped is chosen using swapLocation[]
						double temp = nextPop[(i*2) + 10*a3 + 30][j];
						nextPop[(i*2) + 10*a3 + 30][j] = nextPop[(i*2+1) + 10*a3 + 30][j];
						nextPop[(i*2+1) + 10*a3 + 30][j] = temp;					
					}
				}		
												
			}		
		}
		
		
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
			double counter = pow(integrateSpecies(vPop, i, 90, 180)*0.0548594, 0.5);
			for (int j = 0; j <= SphHarMAX - 1; j++){
				nextPop[i][j] = nextPop[i][j] / counter;
			}		
		}
		
				
// Finally, we equate the old pop to the new pop, allowing the next loop to operate on the new population
		for (int i = 0; i <= PopMAX-1; i++){
			for (int j = 0; j <= SphHarMAX-1; j++){
			pop[i][j] = nextPop[i][j];
			}
		}
					 
	}	// End of Evolution





// The work is done, now time to see the results of the final generation! We follow the same ranking protocol
	double finalScores[100];
	
	for(int i=0; i<PopMAX; i++){
		for(int j=0; j<SphHarMAX; j++){
			vPop[i][j] = pop[i][j];
		}
	}

// Now, to test the curent population, we integrate each species from [theta-spread/2] to [theta + spread/2]
	for (int i = 0; i <= PopMAX - 1; i++){
		finalScores[i] = integrateSpecies(vPop, i, Theta, Spread)*0.0548594;
	}	
	
// We want to organize the population by their testScores, so we make a temporary matrix to hold the ranked initial population:rankedPop 
// Also, a temporary ordered testScores array: rankedTestScores
	double rankedPop[PopMAX][SphHarMAX] = {};
	double rankedFinalScores[PopMAX] = {};
// First, equate the test scores with the ranked test scores. Then we sort the ranked test scores: greatest to least
	for (int i=0; i <= PopMAX-1; i++){
		rankedFinalScores[i] = finalScores[i];
	}

	insertionSort(rankedFinalScores, PopMAX);
		
// Next, we find where testScores = rankedTestScores to place the pop species in the right place in rankedPop
	for (int i=0; i <=PopMAX-1; i++){
		for (int k=0; k <=PopMAX-1; k++){
			if (finalScores[k] == rankedFinalScores[i]){	
				for (int l=0; l <=SphHarMAX-1; l++){ // We need to copy the whole of the species in the kth position to rankedPop in the ith position
					rankedPop[i][l] = pop[k][l];
				}
					
			}
		}
	}

// We print out the highest ranking species's scores and it's arrays in mathematica format, for easy plotting	
	cout << endl << endl << "Final Results:" << endl;

	for (int i = 0; i <= 4; i++){
		cout << "# " << i+1 << " In Mathematica Format with score "<< rankedFinalScores[i]<< " : " << endl; 
		for (int j = 0; j <= SphHarMAX-1; j++){
			cout << "m1Y0" << j << " = " << rankedPop[i][j] << endl;
		}
		cout << endl;
	}	
/*		
// We also print out the whole population, to see the spread of the scores:
	cout << endl << endl << "Final Population:" << endl;

	for (int i = 0; i <= PopMAX-1; i++){
		cout << "Species " << i << " with score "<< rankedFinalScores[i]<< " : " << endl; 
		for (int j = 0; j <= SphHarMAX-1; j++){
			cout << rankedPop[i][j] << "  ";
		}
		cout << endl;
	}		

*/


return 0;
}
