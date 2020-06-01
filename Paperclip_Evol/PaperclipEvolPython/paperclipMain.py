# Written by: Suren Gourapura & Thomas Sinha
# Date: 10/28/17
# Goal:     The goal is to evolve paperclip antennas using rotations as the genetics
#           These antenna will be maximizing curlyness about the z axis

"""
 The following code will evolve paperclip antennas, constructed from unit length segments to maximize different fitness functions
 There are two fitness functions for you to try out for yourself (curliness, and z-height), and one left blank for you to try and write your own
 There will be a set of functions, directly below, that are used to assort different arrays and calculate various fitness scores.
 Below that, in the main, is the code that evolves the antennas. It starts with a random population of 100 different antennas, and then uses
 four different algorithms to maximize a chosen fitness score. More will be explained in the commenting in the main.
"""
import numpy as np
import paperclipFitnessScores
import paperclipGenAlgorithms
import time

def GenPop(popMax, numSeg):
	# np.random.random_sample creates a uniformly random matrix of size (popMax, numSeg, 3) from [0,1)
	population = np.random.random_sample((popMax, numSeg, 3))
	# We want these rotations to be random from the interval [0, 2pi)
	population *= np.pi*2.0
	return population

# The computationally main part of the code. Note: fitType accepts ints from range [1,3]
def RotationMain(numSeg, gen, fitType, fitBreakdown, popMax):
	# Generate the initial population (gen 0)
	pop = GenPop(popMax, numSeg)
	
	# Keep a counter for the number of repeated individuals
	repeatCounter = 0


	for g in range(gen):
		# First, we want to know how well our generation is doing.
		# This function orders them from best to worst and gives a 1D array of scores. 
		rankedScores, rankedPop = paperclipFitnessScores.FitnessTest(pop, fitType)
	
		# Print the best of each generation
		if g % 10==0:
			print 'Generation: ', g
			print "best one",rankedScores[0], rankedPop[0]
			print " "

		"""
		 We first do Algorithm 0.If two antenna have the same fitness score, replace one of them with a random antenna.
		Also, keep a counter of the number of repeated individuals.
		"""
		for i in range(popMax-1):
			if np.abs(rankedScores[i]-rankedScores[i+1])<=0.00001:
				#print "here",i
				repeatCounter += 1
				# Create the new individual that will take the place of the repeated individual
				newIndiv= np.random.random_sample((numSeg, 3))*np.pi*2.0
				# Reshape him into a 1-member population to feed into the ranking code
				newIndivMatrix = newIndiv.reshape((1,numSeg, 3))
				# Feed the 1-member population into the ranking code. The ranked pop here doesn't matter.
				newIndivScore, _ = paperclipFitnessScores.FitnessTest(newIndivMatrix, fitType)
				# Insert the new individual in the population and rank
				rankedPop[i+1] = newIndiv
				rankedScores[i+1] = newIndivScore

		# Re-rank the members to prepare for the remaining algorithms
		rankedScores, rankedPop = paperclipFitnessScores.Sort(rankedScores, rankedPop)	


		# We now give the task of creating the new species to the remaining algorithms individually
		# They will each create the right amount of individuals to sum to popMax
		newPopA1 = paperclipGenAlgorithms.Alg1(rankedScores, rankedPop, fitBreakdown[0])
		newPopA2 = paperclipGenAlgorithms.Alg2(rankedScores, rankedPop, fitBreakdown[1])
		newPopA3 = paperclipGenAlgorithms.Alg3(rankedScores, rankedPop, fitBreakdown[2])
		newPopA4 = paperclipGenAlgorithms.Alg4(rankedScores, rankedPop, fitBreakdown[3])

		newPop = np.vstack((newPopA1, newPopA2, newPopA3, newPopA4))

		pop = newPop
	# End of loop
	

	# Rank the scores one last time
	rankedScores, rankedPop = paperclipFitnessScores.FitnessTest(pop, fitType)
	# Calculate the average repeats per generation
	avgRepeats = (repeatCounter +0.0)/ gen
	# Feed this data to an outputting and formatting function
	OutputFormatting(rankedScores, rankedPop, avgRepeats)

	return

def OutputFormatting(scores, pop, avgRepeats):

	print "Average number of Repeats:", avgRepeats
	
	# Convert the population from rotations to cartesian. Note the extra coordinate introduced 
	cartPop = np.zeros((pop.shape[0], pop.shape[1]+1, pop.shape[2]))
	for i in range(pop.shape[0]):
		cartPop[i] = paperclipFitnessScores.cartTransform(pop[i])
	
	# Now output the data for the best antenna in a Mathematica-friendly format
	print "Antenna with highest fitness score is:"
	print "line = Line [{"
    	for i in range(cartPop.shape[1]):
		if i < cartPop.shape[1]-1:
        		print "{", cartPop[0,i,0] , ", " , cartPop[0,i,1]  ,  ", ", cartPop[0,i,2] ,  "}, "
		if i == cartPop.shape[1]-1:
    			print "{", cartPop[0,i,0] , ", " , cartPop[0,i,1]  ,  ", ", cartPop[0,i,2] , "}}] " 

	print "With a fitness score of", scores[0]


# Run Code Here

# The breakdown controls how many offspring each algorithm is allowed to create. 
# Each value must be divisible by 10 (e.g. 30, not 25)
breakdown = [10, 30, 30, 30]

start = time.time()

RotationMain(numSeg=8, gen=100, fitType=1, fitBreakdown = breakdown, popMax=np.sum(breakdown))

end = time.time()
print ' '
print "The program took:", (end-start), "seconds to run."
print '  '




