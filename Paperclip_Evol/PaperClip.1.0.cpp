
// Written by: Suren Gourapura
// Date: 10/28/17
// Goal: 	The goal is to evolve paperclip antennas using rotations as the genetics
// 			These antenna will be maximizing curlyness about the z axis		

#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
using namespace std;

double Pi = 3.14159265359;
int PopMAX = 100;

// Code modified from https://codereview.stackexchange.com/questions/110793/insertion-sort-in-c
void insertionSort(double array[], int length){ 		// Sort array into greatest -> least
    int i,j;
    for (i = 1; i < length; i++) {
        double temp = array[i];
        for (j = i; j > 0 && array[j - 1] < temp; j--) {
            array[j] = array[j - 1];
        }
        array[j] = temp;
    }
}

void CoordTransform(double oldVec[], double rotx, double roty, double rotz, double newVec[]){
// We are calculating the next unit vector using the previous unit vector and rotations
	double x = oldVec[0];
	double y = oldVec[1];
	double z = oldVec[2];
// The below was calculated on mathematica using the 3D rotation matrices found here: https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
// We are performing rotx first, then roty, then rotz
	newVec[0] = sin(rotx)*(y*sin(roty)*cos(rotz) + z*sin(rotz)) + cos(rotx)*(z*sin(roty)*cos(rotz) - y*sin(rotz)) + x*cos(roty)*cos(rotz);
	newVec[1] = sin(rotz)*(y*sin(rotx)*sin(roty) + x*cos(roty)) + cos(rotx)*(z*sin(roty)*sin(rotz) + y*cos(rotz)) - z*sin(rotx)*cos(rotz);
	newVec[2] = y*sin(rotx)*cos(roty) + z*cos(rotx)*cos(roty) - x*sin(roty);
}

void CrossProduct(double oldVec[], double newVec[], double crossVec[]){
	double x = oldVec[0];
	double y = oldVec[1];
	double z = oldVec[2];
	
	double a = newVec[0];
	double b = newVec[1];
	double c = newVec[2];	
	
	crossVec[0] = c*y - b*z;
	crossVec[1] = a*z - c*x;
	crossVec[2] = b*x - a*y;
}

double FScore(int numSeg, double rotx[], double roty[], double rotz[]){
// First, we need to convert from rotations to unit vector coordinates. 
// Each unit vector corresponds to the direction that line segment is pointing relative to fixed coordinates
// Initialize the converted array
	double unitVecs[numSeg][3];

// The first node is rotated from vertical (0,0,1)
	double newVec[3];
	double oldVec[] = {0,0,1};
	CoordTransform(oldVec, rotx[0], roty[0], rotz[0], newVec);
	for (int i = 0; i < 3; i++){
		unitVecs[0][i] = newVec[i];
	}
	
// Now, we use loops to generate the rest of the coordinates
	for (int i = 1; i < numSeg; i++){
		for (int j = 0; j < 3; j++){									// Get the old vector (the point before point i)
			oldVec[j] = unitVecs[i - 1][j];
		}
		CoordTransform(oldVec, rotx[i], roty[i], rotz[i], newVec);	// Calculate the new vector based on old vector and rotations
		for (int j = 0; j < 3; j++){
			unitVecs[i][j] = newVec[j];								// Add new vectors onto cartesian
		}
	}	
	
// With [numSeg] unit vectors, we can calculate [numSeg-1] cross product vectors.
// The magnitude of these vectors in the z direction (arbitrary choice) gives us the fitness score.
	double crossVec[numSeg - 1][3];
	
	for (int i = 0; i < numSeg - 1; i++){
		for (int j = 0; j < 3; j++){
			oldVec[j] = unitVecs[i][j];		// Initialize the old vector
		}
		for (int j = 0; j < 3; j++){
			newVec[j] = unitVecs[i + 1][j];	// Initialize the new vector
		}	
		CrossProduct(oldVec, newVec, crossVec[i]);
	}

// Simply sum the z component of the cross product vectors to get the fitness score
	double fScore = 0;
	
	for (int i = 0; i < numSeg - 1; i++){
		fScore += crossVec[i][2];
	}
	
	return fScore;
}

double RotToCartesian(int numSeg, double rotx[], double roty[], double rotz[], double xcoord[], double ycoord[], double zcoord[]){
// First, we need to convert from rotations to unit vector coordinates. 
// Each unit vector corresponds to the direction that line segment is pointing relative to fixed coordinates
// Initialize the converted array
	double unitVecs[numSeg][3];

// The first node is rotated from vertical (0,0,1)
	double newVec[3];
	double oldVec[] = {0,0,1};
	CoordTransform(oldVec, rotx[0], roty[0], rotz[0], newVec);
	for (int i = 0; i < 3; i++){
		unitVecs[0][i] = newVec[i];
	}
	
// Now, we use loops to generate the rest of the coordinates
	for (int i = 1; i < numSeg; i++){
		for (int j = 0; j < 3; j++){									// Get the old vector (the point before point i)
			oldVec[j] = unitVecs[i - 1][j];
		}
		CoordTransform(oldVec, rotx[i], roty[i], rotz[i], newVec);	// Calculate the new vector based on old vector and rotations
		for (int j = 0; j < 3; j++){
			unitVecs[i][j] = newVec[j];								// Add new vectors onto cartesian
		}
	}
	
// Now, we need to convert the unit vectors into the actual coordinates. The first point is {0,0,0} and the second point is the same as the first unit vector
// Additional coordinates are made by adding the unit vector onto the previous coordinate
	xcoord[0] = 0;	
	ycoord[0] = 0;
	zcoord[0] = 0;
	
	for (int i = 0; i < numSeg + 1; i++){ // I honestly have no idea why it is numSeg + 1 instead of numSeg, but the latter doesn't work!
		xcoord[i+1] = xcoord[i] + unitVecs[i][0];
		ycoord[i+1] = ycoord[i] + unitVecs[i][1];
		zcoord[i+1] = zcoord[i] + unitVecs[i][2];		
	}	
}




int main()
{
	
	int numSeg;
	int Gen = 0;
	srand(time(NULL));
	
// Enter the number of antenna segments
	cout << "Enter the number of antenna segments, and generations in the following format: #segments generations" << endl;
	cin >> numSeg;
	cin >> Gen;


// Create the population with user specified number of line segments
	double pop[PopMAX][numSeg][3];	
	
// Randomly fill the population. First step is to give each x, y, and z rotation a random value between 0 and 2 pi.
	for (int i = 0; i < PopMAX; i++){       								// for each antenna
		for (int j = 0; j < numSeg; j++){  									// for each segment
			for (int k = 0; k < 3; k++){      								// for each x, y, and z rotation
				pop[i][j][k] = ((double)rand() / (double)(RAND_MAX))*2*Pi; 	// gives a random value between 0 and 2 pi
			}
		}
	}
	

// Now, we evolve. We loop the following code for each generation
	for (int g = 1; g <= Gen; g++){
	//	cout << endl << "Generation: " << g << endl;
		
// The first objective is to calculate how well our current population's species are doing
		double testScores[PopMAX] = {};

// To test the curent population, we first reorganize the 2D array (numSeg x 3) into three 1D arrays, each containing all values for one rotation
		double rotx[numSeg], roty[numSeg], rotz[numSeg];
	
		for (int i = 0; i < PopMAX; i++){       						
			for (int j = 0; j < numSeg; j++){  								
				rotx[j] = pop[i][j][0];
				roty[j] = pop[i][j][1];
				rotz[j] = pop[i][j][2];
			}
			testScores[i] = FScore(numSeg, rotx, roty, rotz); // Each antenna gets its score recorded. The rot arrays are rewritten for each species
		}	
			
// The following matrix will hold the next generation's species. In the end of the evolution, we will make pop = nextPop, so the next evolution afterwards will act on nextPop 
		double nextPop[PopMAX][numSeg][3] = {};



// Evolution Algorithim 1: Take the 10 species with the best scores and pass them onto nextPop

// We want to organize the population by their testScores, so we make a temporary matrix to hold the ranked initial population:rankedPop 
// Also, a temporary ordered testScores array: rankedTestScores
		double rankedPop[PopMAX][numSeg][3] = {};
		double rankedTestScores[PopMAX] = {};
		
// First, equate the test scores with the ranked test scores. Then we sort the ranked test scores: greatest to least
		for (int i = 0; i < PopMAX; i++){
			rankedTestScores[i] = testScores[i];
		}

		insertionSort(rankedTestScores, PopMAX);
		
// Next, we find where testScores = rankedTestScores to place the pop species in the right place in rankedPop
		for (int i = 0; i < PopMAX; i++){
			for (int j = 0; j < PopMAX; j++){
				if (testScores[j] == rankedTestScores[i]){
					
					for (int k = 0; k < numSeg; k++){ // We need to copy the whole of the species in the jth position to rankedPop in the ith position
						for (int l = 0; l < 3; l++){
							rankedPop[i][k][l] = pop[j][k][l];
						}
					}
					
				}
			}
		}
	
// We print out the highest ranking species's score and it's array in mathematica format, for easy plotting	
		cout << rankedTestScores[0] << endl;	
		
// Finally, we copy over the top 10 best species from pop to nextPop		
		for (int i = 0; i < 10; i++){
			for (int j = 0; j < numSeg; j++){
				for (int k = 0; k < 3; k++){
					nextPop[i][j][k] = rankedPop[i][j][k];
				}
			}
		}




		cout << "Generation " << g << " top 5:" << endl;
		for (int i = 0; i < 5; i++){
			cout << "# " << i+1 << " With score "<< rankedTestScores[i]<< " : " << endl; 
			for (int j = 0; j < numSeg; j++){
				cout << "Node " << j << " : X-Rot = " << rankedPop[i][j][0];
				cout << ", Y-Rot = " << rankedPop[i][j][1];
				cout << ", Z-Rot = "<< rankedPop[i][j][2] << endl;	
			}
			cout << endl;
		}
		cout << "# " << 100 << " With score "<< rankedTestScores[99]<< " : " << endl; 
		for (int j = 0; j < numSeg; j++){
			cout << "Node " << j << " : X-Rot = " << rankedPop[99][j][0];
			cout << ", Y-Rot = " << rankedPop[99][j][1];
			cout << ", Z-Rot = "<< rankedPop[99][j][2] << endl;	
		}


// Evolution Algorithim 2: Take 10 random species, find the one with the best score, and [randomly mutate one of it's rotations to obtain an offspring] 10 times. 
// Do this whole algorithim 3 times
		for(int a2 = 1; a2 <= 9; a2++){		
		
			int choose10[10]; // Create an array with 10 random values, 0-99. This array determines the 10 random species that will undergo a tournament selection (a.k.a. simply choosing the highest score species out of the 10)
			for (int i = 0; i < 10; i++){
				choose10[i]= rand() % 100;
			}
		
// Since we have pop already organized from best to worst, we simply find the lowest value in choose10 to find the winner of the tournament. The winning species is called: Algorithim 2 Best Value
			int Alg2BestVal = choose10[0]; // Assume the best species is the first one
		
			for (int i = 1; i < 10; i++){
				if (Alg2BestVal > choose10[i]){  // If another species in the choose10 array has a lower value, it is now the best species
					Alg2BestVal = choose10[i];		
				}
			}
			
// Now, we mutate one part of the species's 2D array, pop[Alg2BestVal][][], to create an offspring. We do this 10 times.
			int mutateLocation[10][2]; // Create the locations for mutation of the 10 offspring. Each location is 2D: which node, which rotation
			for (int i = 0; i < 10; i++){
				mutateLocation[i][0] = rand() % numSeg;	// What segment will be mutated
				mutateLocation[i][1] = rand() % 3;		// What rotation will be mutated
			}
			
// We do the process below 10 times, (to spots 10*a2 to 19*a2 in nextPop)
			for (int i = 0; i < 10; i++){
				for (int j = 0; j < numSeg; j++){
					for (int k = 0; k < 3; k++){
						nextPop[i + 10*a2][j][k] = rankedPop[Alg2BestVal][j][k];	// Copy the Alg2BestVal species over. Note: a2 is added to be able to do the whole of Algorithm 2, 3 times
					}
				}
				for (int j = 0; j < numSeg; j++){
					if (mutateLocation[i][0] == j){			// If the node location is right
						for (int k = 0; k < 3; k++){
							if (mutateLocation[i][1]=k){	// If the rotation location is right
								nextPop[i + 10*a2][j][k] = ((double)rand() / (double)(RAND_MAX))*2*Pi;  // Mutate this location
							}
						}		
					}
				}	
			}			
		}

		
/*		
// Evolution Algorithim 3: Take 20 random species, run two seperate tournaments to find 2 parents, and [swap a random array location with each other to obtain two offspring (for both combinations)] 5 times.
// We do this Algorithm 5 times, for a total of 50 offspring		
		for(int a3 = 1; a3 <= 5; a3++){		
			
			int choose10A[10], choose10B[10]; // Create 2 arrays with 10 random values, 0-99. These arrays determine the 20 random species that will undergo two seperate tournament selections
			int Alg3BestValA = 0, Alg3BestValB = 0;
			
			while (Alg3BestValA == Alg3BestValB){ // To make sure that the 2 parents aren't the same species, we run the following as long as they are the same
				for (int i = 0; i < 10; i++){
					choose10A[i]= rand() % 100;
					choose10B[i]= rand() % 100;
				}
		
// Since we have pop already organized from best to worst, we simply find the lowest value in choose10 to find the winners of the tournaments. The winning species are called: Algorithim 2 Best Value A/B
				Alg3BestValA = choose10A[0]; // Assume the best species are the first ones
				Alg3BestValB = choose10B[0];
				for (int i = 1; i < 10; i++){
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
			for (int i = 0; i < 5; i++){
				swapLocation[i] = rand() % SphHarMAX;
			}
	
// We do the process below 5 times, (to spots 40*a3 to 49*a3 in nextPop)
			for (int i = 0; i < 5; i++){
// First, we copy over the two parents in the 10 offspring spots in nextPop, in an A,B,A,B,... pattern			
				for (int j = 0; j < SphHarMAX; j++){
					nextPop[(i*2) + 10*a3 + 30][j]= rankedPop[Alg3BestValA][j];	// Copy the Alg3BestVal[A and B] species over. Note: a3 is added to be able to do the whole of Algorithm 3, 5 times
					nextPop[(i*2 + 1) + 10*a3 + 30][j]= rankedPop[Alg3BestValB][j];
				}
				for (int j = 0; j < SphHarMAX; j++){
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
		for (int i = 0; i < 10; i++){
			for (int j = 0; j < SphHarMAX; j++){
				nextPop[90 + i][j] = ((double)rand() / (double)(RAND_MAX))*2- 1; // gives a random value between -1 and 1
			}
		}		
	
// Now we normalize nextPop: We divide all coefficients by their integral over all space to normalize to 1
		integral = 0;
		for(int i =0; i < PopMAX; i++){
			for(int j =0; j < SphHarMAX; j++){
				integral += pow(nextPop[i][j], 2);
			}
			for(int j =0; j < SphHarMAX; j++){
				nextPop[i][j] = nextPop[i][j] / sqrt(integral);
			}	
			integral = 0;
		}		
*/				
// Finally, we equate the old pop to the new pop, allowing the next loop to operate on the new population
		for (int i = 0; i < PopMAX; i++){
			for (int j = 0; j < numSeg; j++){
				for (int k = 0; k < 3; k++){
					pop[i][j][k] = nextPop[i][j][k];				
				}
			}
		}
			 
	}	// End of Evolution



// The work is done, now time to see the results of the final generation! We follow the same ranking protocol
	double finalScores[100];


// Now, to test the curent population	
	double rotx[numSeg], roty[numSeg], rotz[numSeg];
	
	for (int i = 0; i < PopMAX; i++){       						
		for (int j = 0; j < numSeg; j++){  								
			rotx[j] = pop[i][j][0];
			roty[j] = pop[i][j][1];
			rotz[j] = pop[i][j][2];
		}
		finalScores[i] = FScore(numSeg, rotx, roty, rotz); // Each antenna gets its score recorded. The rot arrays are rewritten for each species
	}	
	
// We want to organize the population by their testScores, so we make a temporary matrix to hold the ranked initial population:rankedPop 
// Also, a temporary ordered testScores array: rankedTestScores
	double rankedPop[PopMAX][numSeg][3] = {};
	double rankedFinalScores[PopMAX] = {};
// First, equate the test scores with the ranked test scores. Then we sort the ranked test scores: greatest to least
	for (int i=0; i < PopMAX; i++){
		rankedFinalScores[i] = finalScores[i];
	}

	insertionSort(rankedFinalScores, PopMAX);
		
// Next, we find where testScores = rankedTestScores to place the pop species in the right place in rankedPop
		for (int i = 0; i < PopMAX; i++){
			for (int j = 0; j < PopMAX; j++){
				if (finalScores[j] == rankedFinalScores[i]){
					
					for (int k = 0; k < numSeg; k++){ // We need to copy the whole of the species in the jth position to rankedPop in the ith position
						for (int l = 0; l < 3; l++){
							rankedPop[i][k][l] = pop[j][k][l];
						}
					}
					
				}
			}
		}
/*
	for (int i = 0; i < PopMAX; i++){
		for (int k = 0; k < PopMAX; k++){
			if (finalScores[k] == rankedFinalScores[i]){	
				for (int l = 0; l < SphHarMAX; l++){ // We need to copy the whole of the species in the kth position to rankedPop in the ith position
					rankedPop[i][l] = pop[k][l];
				}
					
			}
		}
	} */

// We print out the highest ranking species's scores and it's arrays in mathematica format, for easy plotting	
	cout << endl << endl << "Final Results:" << endl;

	for (int i = 0; i < 5; i++){
		cout << "# " << i+1 << " With score "<< rankedFinalScores[i]<< " : " << endl; 
		for (int j = 0; j < numSeg; j++){
			cout << "Node " << j << " : X-Rot = " << rankedPop[i][j][0];
			cout << ", Y-Rot = " << rankedPop[i][j][1];
			cout << ", Z-Rot = "<< rankedPop[i][j][2] << endl;	
		}
		cout << endl;
	}
	
// To print in Mathematica friendly format, we first need to convert the rotations into unit vectors, then into cartesian coordinates	



	cout << endl << endl << "Final Results in Mathematica Format:" << endl;

	for (int i = 0; i < 5; i++){
		cout << "# " << i+1 << " With score "<< rankedFinalScores[i]<< " : " << endl;
		double rotx[numSeg], roty[numSeg], rotz[numSeg];
		double xcoord[numSeg+1], ycoord[numSeg+1], zcoord[numSeg+1];
							
		for (int j = 0; j < numSeg; j++){  								
			rotx[j] = rankedPop[i][j][0];
			roty[j] = rankedPop[i][j][1];
			rotz[j] = rankedPop[i][j][2];
		}		
	
		RotToCartesian(numSeg, rotx, roty, rotz, xcoord, ycoord, zcoord);

		cout << "line = Line [{";
   		for (int i = 0; i < numSeg; i++){
  	 		cout << "{" << xcoord[i] << ", " << ycoord[i] << ", "<< zcoord[i] << "}, ";
		}
		cout << "{" << xcoord[numSeg] << ", " << ycoord[numSeg] << ", "<< zcoord[numSeg] << "}}] " << endl;
		}
return 0;
}
