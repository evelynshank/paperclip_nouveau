

import numpy as np

def Tournament(rPop, numCompetitors):
	# Gives (the index of) 1 winner out of a tournament of [numCompetitors] number of competitors
	# Create an array with all indicies in rPop: [0,1,2,..., popMax]
	totalPopInd = np.arange(rPop.shape[0])
	# Shuffle it randomly
	np.random.shuffle(totalPopInd)
	# Choose the first numCompetitors
	tournCompetitors = totalPopInd[:numCompetitors]
	# The competitor with the best fitness score is the competitor with the lowest number, since the population is ordered
	return np.amin(tournCompetitors)






# Survival of the Fittest
def Alg1(rScores, rPop, numOffspring):
	"""
	The easiest algorithm! Send back the top numOffspring number of the best individuals
	"""
	return rPop[:numOffspring]

# Mutation (Asexual Reproduction)
def Alg2(rScores, rPop, numOffspring):
	"""
	We take 10 random species, find the one with the best score, and randomly mutate one of its rotations 10 seperate times to obtain 10 offspring. This whole process is done (numOffspring/10) times.
	"""
	offspring = np.zeros((numOffspring, rPop.shape[1], rPop.shape[2]))
	
	for i in range(numOffspring/10):
		# We need to perform a tournament selection first to get 1 winner out of 10
		parent = rPop[Tournament(rPop, 10)]
	
		for j in range(10):
			# We need a location for mutation (node and x,y,or z) and a mutation value
			whichNode = np.random.randint(rPop.shape[1])
			whichRot = np.random.randint(3)
			newVal = np.random.random_sample()*np.pi*2.0
			# Now we put these values into the new offspring
			offspring[i*10+j] = parent
			offspring[i*10+j, whichNode, whichRot] = newVal
	
	return offspring





# Fine Mutation
def Alg3(rScores, rPop, numOffspring):
	"""
	We take 10 random species, find the one with the best score, and randomly mutate one of its rotations 10 seperate times to obtain 10 offspring. This time, we simply tune a rotation up or down slightly instead of replacing it entirely. This whole process is done (numOffspring/10) times.
	"""
	offspring = np.zeros((numOffspring, rPop.shape[1], rPop.shape[2]))
	
	for i in range(int(numOffspring/10.0)):
		# We need to perform a tournament selection first, to get 1 winner out of 10
		parent = rPop[Tournament(rPop, 10)]

		for j in range(10):
			# We need a location for mutation (node and x,y,or z)
			whichNode = np.random.randint(rPop.shape[1])
			whichRot = np.random.randint(3)
			# This time, we also need a fine tuning factor epsilon
			epsilon = np.pi/6.0
			newVal = parent[whichNode, whichRot] + epsilon
			# We may run into a problem here, when newVal > 2 pi. So, adjust for this
			if newVal >= 2*np.pi:
				newVal -= 2*np.pi
			# Now we put these values into the new offspring
			offspring[i*10+j] = parent
			offspring[i*10+j, whichNode, whichRot] = newVal
	return offspring

# Crossover (Sexual Reproduction)
def Alg4(rScores, rPop, numOffspring):
	"""
	We take two sets of 10 random species, find the best score in each, and randomly crossover one of its rotations five seperate times to obtain 10 offspring. Each switch between parents A and B creates two offspringone made of mostly A and one made with mostly B. This whole This whole process is done (numOffspring/10) times. 
	"""
	offspring = np.zeros((numOffspring, rPop.shape[1], rPop.shape[2]))
	
	for i in range(int(numOffspring/10.0)):
		# We need to perform two tournament selections first. Gets 2 parents that aren't the same
		parentAind = Tournament(rPop, 10)
		parentBind = Tournament(rPop, 10)
		while parentAind == parentBind:
			#print 'repeated parent', i
			parentBind = Tournament(rPop, 10)
		parentA = rPop[parentAind]
		parentB = rPop[parentBind]
		for j in range(5):
			# We need a location for crossover (node and x,y,or z)
			whichNode = np.random.randint(rPop.shape[1])
			whichRot = np.random.randint(3)
			# This process creates two offspring. First, copy the parents over
			offspring[i*10+j*2] = parentA
			offspring[i*10+j*2+1] = parentB
			# Now, switch the specific rotation in each offspring
			offspring[i*10+j*2, whichNode, whichRot] = parentB[whichNode, whichRot]
			offspring[i*10+j*2+1, whichNode, whichRot] = parentA[whichNode, whichRot]
	return offspring











