"""
   Multifractal Deterended Fluctuation Analysis (MFDFA)
   REF : Kantelhardt et al., Physica A, 316, 1, 2002.
   
   Original Octave code is translated in Python by Neelakshi Joshi.
   
   Few steps have been adapted from author "Dominik Krzeminski (dokato)" as it gives easier way to cut the series in smaller segments to do    
   local detrending. Auhtor's original code for DFA can be found at https://github.com/dokato/dfa/blob/master/dfa.py
   
   Inputs: 
	data = timeseries (1-d array)
	m = detrending order (int)
        scale = scales (1-d array)
        q = moments (1-d array)
             needCorrection parameter decides how to define profile to implement the MFDFA, 
             according to the conditions described by Ihlen, 2012 in table 1.
	needCorrection = 0 then profile = cumsum(data - mean(data)
	needCorrection = 1 then profile = diff(data)	
	needCorrection = value other than 0 and 1 then profile = data
   
   Outputs:
	Fqs = qth order Fluctuation function
	hq = the generalized Hurst exponents
	tq = the classical multifractal exponents
	alpha and falpha = multifractal spectrum

"""

import numpy as np

def mfdfa_frwd(data, order, scales, q, needCorrection):

    N = len(data)
    if needCorrection == 0:
        profile = np.cumsum(data - np.mean(data))
    if needCorrection == 1:
        profile = np.diff(data)
   
    Fqs =  np.zeros([len(q),len(scales)])
    
    for index_sc, sc in enumerate(scales):
        seg = (data.shape[0] // sc, sc)
        X = np.lib.stride_tricks.as_strided(data, shape = seg)
        scale_ax = np.arange(sc)
        rms = np.zeros(X.shape[0])
        for indx, series_seg in enumerate(X):
            coeff = np.polyfit(scale_ax, series_seg, 1)
            xfit = np.polyval(coeff, scale_ax)
            rms[indx] = np.sqrt(np.mean((series_seg - xfit) ** 2))

        for index_q, iq in enumerate(q):
            if iq == 0:
                Fqs[index_q,index_sc] = np.exp(0.5 * np.mean(np.log(rms**2)))
            else:
                qrms = rms ** iq
                Fqs[index_q, index_sc] = np.mean(qrms) ** (1 / iq)

       
    hq = np.zeros(len(q))
    tauq = np.zeros(len(q))
    alpha = np.zeros(len(q)-1)
    falpha = np.zeros(len(q)-1)

    for num_q in np.arange(len(q)):
        coeff = np.polyfit(np.log2(scales), np.log2(Fqs[num_q,]), 1)
        hq[num_q] = coeff[0]
        
    tauq = hq * q - 1
    alpha = np.diff(tauq) / (q[1]-q[0])
    falpha = (q[:-1]*alpha) - tauq[:-1]
    
    return Fqs, hq, tq, alpha, falpha


if __name__ == "__main__":
    
     Fqs, hq, tq, alpha, falpha = mfdfa_frwd(data, order, scales, q, needCorrection)
