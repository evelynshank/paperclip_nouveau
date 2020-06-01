This code was written in late 2017.

The goal was to write code that would model a gaussian wrapped around a point. This is to simulate a gain pattern that is sensitive to a single direction

The inputs are: theta phi stdev

The code develops a noprmalized gaussian with user inputted standard deviation centered on user inputted theta and phi. The output is a matrix for the height of the gaussian/gain function at each theta and phi in 5 degree increments.

Provided is a Mathematica notebook to copy and paste the matrix and plot the result against a normalized sphere.