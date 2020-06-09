# codes
This repository contains codes for MFDFA and P model in Octave as well as in Python. 
These codes are part of my thesis work titled 
"Multifractal and cascade P model analysis for characterizing non-homogeneous turbulence in space physics" 
done at National Institute for Space Research (INPE), Brazil; 
under the supervision of Dr. Reinaldo Rosa & Dr. Stephan Stephany.

mfdfa1d: mutifractal analysis method viz Multifractal Detrended Fluctuation Analysis (MFDFA) is proposed by 
Kantelhardt et al. (2002) to characterize multiple scaling behavior in the data.
MFDFA is implemented in Octave and in Python. The code works on one-dimensional time series.

pmodel: P model is proposed by Meneveau and Sreenivasan (1987) which is based on the generalized two–scale Cantor set 
to mimic the possible energy transfer in the turbulent cascade process. The P model always describes cascading processes in the
multifractal time series. The P model considers equal lengths (l1 = l2) and unequal weights (p1 != p2 and p1 + p2 +dp = 1). 
For further details refer: 
Joshi et al., Structural characterization of the equatorial F region plasma irregularities in the multifractal context,
Annales Geophysicae, 38, 2, 2020. (https://www.ann-geophys.net/38/445/2020/)
Please cite the above mentioned article when using this code. 
