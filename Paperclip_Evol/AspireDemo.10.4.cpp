// Written by: Suren Gourapura & Thomas Sinha
// Date: 10/28/17
// Goal:     The goal is to evolve paperclip antennas using rotations as the genetics
//           These antenna will be maximizing curlyness about the z axis

/*
 The following code will evolve paperclip antennas, constructed from unit length segments to maximize different fitness functions
 There are two fitness functions for you to try out for yourself (curliness, and z-height), and one left blank for you to try and write your own
 There will be a set of functions, directly below, that are used to assort different arrays and calculate various fitness scores.
 Below that, in the main, is the code that evolves the antennas. It starts with a random population of 100 different antennas, and then uses
 four different algorithms to maximize a chosen fitness score. More will be explained in the commenting in the main.
 */

#include <iomanip>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
#include <fstream>
using namespace std;

double Pi = 3.14159265359;
const int PopMAX = 100; //100
int Fit;


// The following are functions used in the calculation of fitness scores. The populations are constructed from rotations, which we can transfer to cartesian coordinates and use Mathematica to plot

// Code modified from https://codereview.stackexchange.com/questions/110793/insertion-sort-in-c
void insertionSort(double array[], int length){         // Sort array into greatest -> least
    int i,j;
    for (i = 1; i < length; i++) {
        double temp = array[i];
        for (j = i; j > 0 && array[j - 1] < temp; j--) {
            array[j] = array[j - 1];
        }
        array[j] = temp;
    }
}

// The classes classes below, namely CoordTransform, CrossProduct, and RotToCartesian are simply mathematical objects used to calculate
// the unit vectors, and the curliness. You will not have to worry about them.

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
        for (int j = 0; j < 3; j++){                                    // Get the old vector (the point before point i)
            oldVec[j] = unitVecs[i - 1][j];
        }
        CoordTransform(oldVec, rotx[i], roty[i], rotz[i], newVec);    // Calculate the new vector based on old vector and rotations
        for (int j = 0; j < 3; j++){
            unitVecs[i][j] = newVec[j];                                // Add new vectors onto cartesian
        }
    }
    
    // Now, we need to convert the unit vectors into the actual coordinates. The first point is {0,0,0} and the second point is the same as the first unit         vector
    // Additional coordinates are made by adding the unit vector onto the previous coordinate
    xcoord[0] = 0;
    ycoord[0] = 0;
    zcoord[0] = 0;
    
    for (int i = 0; i < numSeg + 1; i++){ // I honestly have no idea why it is numSeg + 1 instead of numSeg, but the latter doesn't work!
        xcoord[i+1] = xcoord[i] + unitVecs[i][0];
        ycoord[i+1] = ycoord[i] + unitVecs[i][1];
        zcoord[i+1] = zcoord[i] + unitVecs[i][2];
    }
    return 0;
}

// This class calculates a fitness score that corresponds mathematically to the discrete curl of the antenna
// All this means is that we take the z-component of the cross-product of every two adjacent vectors, and
// sum them up for the fitness score. An ideal antenna will be coplanar with the x-y plane, and will be
// a square coil
double FScore1(int numSeg, double rotx[], double roty[], double rotz[]){
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
        for (int j = 0; j < 3; j++){                                    // Get the old vector (the point before point i)
            oldVec[j] = unitVecs[i - 1][j];
        }
        CoordTransform(oldVec, rotx[i], roty[i], rotz[i], newVec);    // Calculate the new vector based on old vector and rotations
        for (int j = 0; j < 3; j++){
            unitVecs[i][j] = newVec[j];                                // Add new vectors onto cartesian
        }
    }
    
    // With [numSeg] unit vectors, we can calculate [numSeg-1] cross product vectors.
    // The magnitude of these vectors in the z direction (arbitrary choice) gives us the fitness score.
    double crossVec[numSeg - 1][3];
    
    for (int i = 0; i < numSeg - 1; i++){
        for (int j = 0; j < 3; j++){
            oldVec[j] = unitVecs[i][j];        // Initialize the old vector
        }
        for (int j = 0; j < 3; j++){
            newVec[j] = unitVecs[i + 1][j];    // Initialize the new vector
        }
        CrossProduct(oldVec, newVec, crossVec[i]);
    }
    
    // Simply sum the z component of the cross product vectors to get the fitness score
    double fScore1 = 0;
    
    for (int i = 0; i < numSeg - 1; i++){
        fScore1 += crossVec[i][2];
    }
    
    return fScore1;
}


// This class calculates a fitness score that measures the z-height of an antenna. It’s very simple and should be the first one you look at
// The ideal antenna from this class would just be a vertical line.

double FScore2 (int numSeg, double rotx[], double roty[], double rotz[]){
    double unitVecs[numSeg][3];
    double newVec[3];
    double oldVec[] = {0,0,1};
    CoordTransform(oldVec, rotx[0], roty[0], rotz[0], newVec);
    for (int i = 0; i < 3; i++){
        unitVecs[0][i] = newVec[i];
    }
    
    for (int i = 1; i < numSeg; i++){
        for (int j = 0; j < 3; j++){                                    // Get the old vector (the point before point i)
            oldVec[j] = unitVecs[i - 1][j];
        }
        CoordTransform(oldVec, rotx[i], roty[i], rotz[i], newVec);    // Calculate the new vector based on old vector and rotations
        for (int j = 0; j < 3; j++){
            unitVecs[i][j] = newVec[j];                                // Add new vectors onto cartesian
        }
    }
    
    double fScore2 = 0;
    for(int i = 0; i < numSeg; i++){
        fScore2 += unitVecs[i][2];//-unitVecs[i][0]-unitVecs[i][1];
    }
    return fScore2;
}
/* And lastly, this class is yours to play with. See what kind of fitness scores you can come up with, and try to implement them.
 Some ideas: Minimize/Maximize length in any of the directions, make the antenna loop back on itself (This one is fun. Think distance formula)
 */
double FScore3 (int numSeg, double rotx[], double roty[], double rotz[]){
    double xcoord[numSeg+1], ycoord[numSeg+1], zcoord[numSeg+1];
    double fScore3 = 0;
    RotToCartesian(numSeg,rotx,roty,rotz,xcoord,ycoord,zcoord);
    fScore3 = 1/(abs(xcoord[numSeg+1])+abs(ycoord[numSeg+1])+abs(zcoord[numSeg+1])+0.01);
    return fScore3;
}
int main()
{
    const int numSeg = 10;
    int Gen = 0;
    srand((unsigned int)time(NULL));
    
    // Enter the number of antenna segments
    cout << "Enter the number of generations:" << endl;
    // cin >> numSeg;
    cin >> Gen;
    cout << endl << "Choose your Fitness Function:" << endl << "1.) Curliness" << endl << "2.) Z-Height" << endl << "3.) Loopiness"<< endl;
    cin >> Fit;
    
    // Create the population with user specified number of line segments
    double pop[PopMAX][numSeg][3];
    
    // Randomly fill the population. First step is to give each x, y, and z rotation a random value between 0 and 2 pi.
    for (int i = 0; i < PopMAX; i++){                                       // for each antenna
        for (int j = 0; j < numSeg; j++){                                      // for each segment
            for (int k = 0; k < 3; k++){                                      // for each x, y, and z rotation
                pop[i][j][k] = ((double)rand() / (double)(RAND_MAX))*2*Pi;     // gives a random value between 0 and 2 pi
            }
        }
    }
    
    
    // Now, we evolve. We loop the following code for each generation
    for (int g = 1; g <= Gen; g++){
        //    cout << endl << "Generation: " << g << endl;
        
        // The first objective is to calculate how well our current population's species are doing
        double testScores[PopMAX] = {};
        
        
        // Here, we can choose which fitness function we wish to use, using an if statement
        
        if (Fit == 1){   // The block of code in this if statement calculates the Curliness fitness score of our population
            // To test the curent population, we first reorganize the 2D array (numSeg x 3) into three 1D arrays, each containing all values for one rotation
            double rotx[numSeg], roty[numSeg], rotz[numSeg];
            
            for (int i = 0; i < PopMAX; i++){
                for (int j = 0; j < numSeg; j++){
                    rotx[j] = pop[i][j][0];
                    roty[j] = pop[i][j][1];
                    rotz[j] = pop[i][j][2];
                }
                testScores[i] = FScore1(numSeg, rotx, roty, rotz); // Each antenna gets its score recorded. The rot arrays are rewritten for each species
            }
        }
        
        else if (Fit == 2){
            for(int i=0; i<PopMAX; i++){
                double rotx[numSeg], roty[numSeg], rotz[numSeg];
                
                for (int i = 0; i < PopMAX; i++){
                    for (int j = 0; j < numSeg; j++){
                        rotx[j] = pop[i][j][0];
                        roty[j] = pop[i][j][1];
                        rotz[j] = pop[i][j][2];
                    }
                }
                
                testScores[i] = FScore2(numSeg, rotx ,roty , rotz);
            }
        }
        else if (Fit == 3){
            for (int i= 0; i< PopMAX; i++){
                double rotx[numSeg], roty[numSeg], rotz[numSeg];
                
                for (int i = 0; i < PopMAX; i++){
                    for (int j = 0; j < numSeg; j++){
                        rotx[j] = pop[i][j][0];
                        roty[j] = pop[i][j][1];
                        rotz[j] = pop[i][j][2];
                    }
                }
                
                testScores[i] = FScore3(numSeg, rotx ,roty , rotz);
            }
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
        
        // Finally, we copy over the top 10 best species from pop to nextPop
        for (int i = 0; i < 10; i++){
            for (int j = 0; j < numSeg; j++){
                for (int k = 0; k < 3; k++){
                    nextPop[i][j][k] = rankedPop[i][j][k];
                }
            }
        }
        
        
        // Evolution Algorithim 2: Take 10 random species, find the one with the best score, and [randomly mutate one of it's rotations to obtain an offspring] 10 times.
        // Do this whole algorithim 4 times
        for(int a2 = 1; a2 <= 5; a2++){
            
            int choose10[10]; // Create an array with 10 random values, 0-99. This array determines the 10 random species that will undergo a tournament selection (a.k.a. simply choosing the highest score species out of the 10)
            for (int i = 0; i < 10; i++){
                choose10[i]= rand() % PopMAX;
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
                mutateLocation[i][0] = rand() % numSeg;    // What segment will be mutated
                mutateLocation[i][1] = rand() % 3;        // What rotation will be mutated
            }
            
            // We do the process below 10 times, (to spots 10*a2 to 19*a2 in nextPop)
            for (int i = 0; i < 10; i++){
                for (int j = 0; j < numSeg; j++){
                    for (int k = 0; k < 3; k++){
                        nextPop[i + 10*a2][j][k] = rankedPop[Alg2BestVal][j][k];    // Copy the Alg2BestVal species over. Note: a2 is added to be able to do the whole of Algorithm 2, 3 times
                    }
                }
                for (int j = 0; j < numSeg; j++){
                    if (mutateLocation[i][0] == j){            // If the node location is right
                        for (int k = 0; k < 3; k++){
                            if (mutateLocation[i][1]==k){    // If the rotation location is right
                                nextPop[i + 10*a2][j][k] = ((double)rand() / (double)(RAND_MAX))*2*Pi;  // Mutate this location
                            }
                        }
                    }
                }
            }
        }
        
        
        
        // Evolution Algorithim 3: Take 20 random species, run two seperate tournaments to find 2 parents, and [swap a random array location with each other to obtain two offspring (for both combinations)] 5 times.
        // We do this Algorithm 5 times, for a total of 50 offspring
        for(int a3 = 1; a3 <= 3; a3++){
            
            int choose10A[10], choose10B[10]; // Create 2 arrays with 10 random values, 0-99. These arrays determine the 20 random species that will undergo two seperate tournament selections
            int Alg3BestValA = 0, Alg3BestValB = 0;
            
            while (Alg3BestValA == Alg3BestValB){ // To make sure that the 2 parents aren't the same species, we run the following as long as they are the same
                for (int i = 0; i < 10; i++){
                    choose10A[i]= rand() % PopMAX;
                    choose10B[i]= rand() % PopMAX;
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
            int swapLocation[numSeg]; // Create the locations for swapping
            for (int i = 0; i < numSeg; i++){
                swapLocation[i] = rand() %  numSeg;
            }
            
            // We do the process below 5 times, (to spots 50*a3 to 59*a3 in nextPop)
            for (int i = 0; i < 5; i++){
                // First, we copy over the two parents in the 10 offspring spots in nextPop, in an A,B,A,B,... pattern
                for (int j = 0; j < numSeg ; j++){
                    for (int k = 0; k<3; k++){
                        nextPop[(i*2) + 10*a3 + 60][j][k]= rankedPop[Alg3BestValA][j][k];    // Copy the Alg3BestVal[A and B] species over. Note: a3 is added to be able to do the whole of Algorithm 3, 5 times
                        nextPop[(i*2 + 1) + 10*a3 + 60][j][k]= rankedPop[Alg3BestValB][j][k];
                    }
                }
                for (int j = 0; j < numSeg ; j++){
                    for (int k = 0; k<3; k++){
                        if(swapLocation[i] == j){  // make sure that the location to be swapped is chosen using swapLocation[]
                            double temp[5]={};
                            temp[k] = nextPop[(i*2) + 10*a3 + 30][j][k];
                            nextPop[(i*2) + 10*a3 + 60][j][k] = nextPop[(i*2+1) + 10*a3 + 30][j][k];
                            nextPop[(i*2+1) + 10*a3 + 60][j][k] = temp[k];
                        }
                    }
                }
                
            }
        }
        
        
        
        // Evolution Algorithim 4: Introduce 10 random species into the population
        // First step is to give each value in each species random number between -1 and 1. Then, we normalize the species values to 1
        for (int i = 0; i < 10; i++){
            for (int j = 0; j < numSeg; j++){
                for(int k = 0; k<3;k++){
                    nextPop[90 + i][j][k] = ((double)rand() / (double)(RAND_MAX))*2*Pi; // gives a random value between -1 and 1
                }
            }
        }
        
        
        // Finally, we equate the old pop to the new pop, allowing the next loop to operate on the new population
        for (int i = 0; i < PopMAX; i++){
            for (int j = 0; j < numSeg; j++){
                for (int k = 0; k < 3; k++){
                    pop[i][j][k] = nextPop[i][j][k];
                }
            }
        }
    }    // End of Evolution
    
    // The work is done, now time to see the results of the final generation! We follow the same ranking protocol
    double finalScores1[PopMAX];
    double finalScores2[PopMAX];
    double finalScores3[PopMAX];
    
    // Now, to test the curent population
    double rotx[numSeg], roty[numSeg], rotz[numSeg];
    
    for (int i = 0; i < PopMAX; i++){
        for (int j = 0; j < numSeg; j++){
            rotx[j] = pop[i][j][0];
            roty[j] = pop[i][j][1];
            rotz[j] = pop[i][j][2];
        }
        finalScores1[i] = FScore1(numSeg, rotx, roty, rotz); // Each antenna gets its score recorded. The rot arrays are rewritten for each species
        finalScores2[i] = FScore2(numSeg, rotx, roty, rotz);
        finalScores3[i] = FScore3(numSeg, rotx, roty, rotz);
    }
    
    // We want to organize the population by their testScores, so we make a temporary matrix to hold the ranked initial population:rankedPop
    // Also, a temporary ordered testScores array: rankedTestScores
    double rankedPop[PopMAX][numSeg][3] = {};
    double rankedFinalScores[PopMAX] = {};
    
    // First, equate the test scores with the ranked test scores. Then we sort the ranked test scores: greatest to least
    for (int i=0; i < PopMAX; i++){
        if (Fit== 1){
            rankedFinalScores[i] = finalScores1[i];
        }
        else if (Fit == 2){
            rankedFinalScores[i] = finalScores2[i];
        }
        else if (Fit == 3){
            rankedFinalScores[i] = finalScores3[i];
        }
    }
    
    insertionSort(rankedFinalScores, PopMAX);
    
    // Next, we find where testScores = rankedTestScores to place the pop species in the right place in rankedPop
    if (Fit == 1){
        for (int i = 0; i < PopMAX; i++){
            for (int j = 0; j < PopMAX; j++){
                if (finalScores1[j] == rankedFinalScores[i]){
                    for (int k = 0; k < numSeg; k++){ // We need to copy the whole of the species in the jth position to rankedPop in the ith position
                        for (int l = 0; l < 3; l++){
                            rankedPop[i][k][l] = pop[j][k][l];
                        }
                    }
                    
                }
            }
        }
    }
    else if (Fit == 2){
        for (int i = 0; i < PopMAX; i++){
            for (int j = 0; j < PopMAX; j++){
                if (finalScores2[j] == rankedFinalScores[i]){
                    for (int k = 0; k < numSeg; k++){ // We need to copy the whole of the species in the jth position to rankedPop in the ith position
                        for (int l = 0; l < 3; l++){
                            rankedPop[i][k][l] = pop[j][k][l];
                        }
                    }
                    
                }
            }
        }
    }
    else if (Fit == 3){
        for (int i = 0; i < PopMAX; i++){
            for (int j = 0; j < PopMAX; j++){
                if (finalScores3[j] == rankedFinalScores[i]){
                    for (int k = 0; k < numSeg; k++){ // We need to copy the whole of the species in the jth position to rankedPop in the ith position
                        for (int l = 0; l < 3; l++){
                            rankedPop[i][k][l] = pop[j][k][l];
                        }
                    }
                    
                }
            }
        }
    }
    
    // To print in Mathematica friendly format, we first need to convert the rotations into unit vectors, then into cartesian coordinates
    
    cout << endl << endl << "Final Results in Mathematica Format:" << endl;
    
    cout << "#1 With score "<< rankedFinalScores[0]<< " : " << endl;
    // double rotx[numSeg], roty[numSeg], rotz[numSeg];
    double xcoord[numSeg+1], ycoord[numSeg+1], zcoord[numSeg+1];
    
    for (int j = 0; j < numSeg; j++){
        rotx[j] = rankedPop[0][j][0];
        roty[j] = rankedPop[0][j][1];
        rotz[j] = rankedPop[0][j][2];
        
        RotToCartesian(numSeg, rotx, roty, rotz, xcoord, ycoord, zcoord);
    }
    cout << "line = Line [{";
    for (int i = 0; i < numSeg; i++){
        cout << "{" << xcoord[i] << ", " << ycoord[i] << ", "<< zcoord[i] << "}, ";
    }
    cout << "{" << xcoord[numSeg] << ", " << ycoord[numSeg] << ", "<< zcoord[numSeg] << "}}] " << endl;
    
    
    return 0;
}
