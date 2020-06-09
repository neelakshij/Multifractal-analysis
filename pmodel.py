#################################################################################################################################
'''
   Credit: Neelakshi Joshi 
   This code is written by Dr. Neelakshi Joshi for her thesis work done at National Institute for Space Research (INPE), Brazil.
   Original Octave code is translated in Python by Neelakshi Joshi.

   This code generates anaytical "P model singularity spectrum" based on the generalized two-scale weighted Cantor set.
   References: 1) Meneveau and Sreenivasan, Physical Review Letters, 59, 1987.
	       2) Halsey et al., Physical Review A, 33, 2, 1986.

   Inputs: 
	p1 = probability parameter, 0 < p1 < 1
	l1 = cascading length parameter, 0 < l1 < 1
	dp = loss of probability where p1+p2 < 1 and p1+p2+dp = 1
   Output:
	alpha and falpha = singularity spectrum

'''
###################################################################################################################################

import numpy as np

def pmodel(p1, l1, dp):
    p2 = (1-p1) - dp
    l2 = 1-l1
    q = np.arange(-20, 20, 0.5)
    alpha = np.zeros(len(q))
    falpha = np.zeros(len(q))

    weight = p2 / p1
    nbym = 1 + (weight) ** q
    mbyn = 1. / nbym
    bt = (weight) ** q

    for i in np.arange(1,len(bt)):
        alpha[i] = (np.log(p1) + (bt[i] * np.log(p2))) / (np.log(l1) + bt[i] * np.log(l2))
        falpha[i] = (bt[i] * np.log(bt[i]) - (bt[i] + 1) * np.log(bt[i] + 1)) / (np.log(l1) + bt[i] * np.log(l2))
    
    return alpha, falpha
 
if __name__ == "__main__":
    
    alpha, falpha = pmodel(0.4, 0.5, 0.0)
