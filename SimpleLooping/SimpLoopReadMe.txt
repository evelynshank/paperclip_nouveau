This code was written on spring break 2018.

The goal was to see how a simple looping code compares to the ring evolution code. 

The code uses the same fitness score from the ring evolution code, but instead of evolving populations of coefficients, it generates 13 random coefficients of the spherical harmonics, checks their fitness value, and if their score is worse than the previous best score, it trys again. If their score is better than the previous best score, the coefficients and fitness score are saved as the new best score and species. Once a species is generated with a fitness score better than what is specified in the code (e.g. a double number from 0 to 1, usually between 0.1 and 0.8), then the code stops. Dev-C++ provides a timestamp for how long the code ran for, so this value can be recoded to judge the effectiveness of simple looping.

Plotting can be done with the same notebook as is used with the ring evolution code. Included is an excel with tabulated results for times comparing the ring evolver with the simple looping code. Note that these trials were run on a laptop, and all times would reduce drastically on a faster machine (factor of 2 or more on my desktop).