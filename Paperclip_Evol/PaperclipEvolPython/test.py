import numpy as np
import paperclipFitnessScores
import paperclipGenAlgorithms

#scores = np.array([4.5,3,6])
#indivs = np.array([[[1,1,1],[2,2,2],[1,1,1],[2,2,2]], [[3,3,3],[4,4,4],[3,3,3],[4,4,4]],[ [5,5,5],[6,6,6],[5,5,5],[6,6,6]]])

vecMatrixIndiv = np.array([ [0,1,2],[5,-4,2],[0,-1,-2],[-5,4,-2] ])
coordMatrixIndiv = np.zeros((vecMatrixIndiv.shape[0]+1,vecMatrixIndiv.shape[1]))

for i in range(coordMatrixIndiv.shape[0]): # <- this should have the value of numSeg
		coordMatrixIndiv[i]=np.sum(vecMatrixIndiv[0:i], axis=0)




print vecMatrixIndiv 

print coordMatrixIndiv 
