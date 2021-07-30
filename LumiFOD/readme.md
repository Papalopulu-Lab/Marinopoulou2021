The LumiFOD folder contains routines developed by Nick Phillips which use a periodic and aperiodic Gaussian process covariance model to infer parameters by
maximizing likelihood and select statistically signficant periodic timeseries based on the log likelihood ratio (LLR) statistic. This version contains minor 
customisations including: optimisation of covariance models separately before computing LLR, routines to combine and export single cell statistics and detrended data.
In addition, the Hilbert transform was used to measure the peak to through fold change in amplitude. 

To use this please cite the following publications:

Phillips et al. (2017) "Identifying stochastic oscillations in single-cell live imaging time series using Gaussian processes" PLoS Comput Biol 13(5): e1005479. 
Here you can find the the covariance models and the overall approach.

C Rasmussen and N Hannes (2010) "Gaussian processes for machine learning (GPML) toolbox" Journal of Machine Learning Research, 11:3011 3015. 
Here you can find the information about the GPML toolbox. The copy shared with this repository contains additional routines. 

To run the code please follow these steps:

- download the GMPL. zip provided in the same repository and unarchive this
- download a copy of the LumiFOD folder in the same location as the above
- run FODmain_exportHilbert.m
- results will be stored in a simulation folder as figures, .mat and .xls files
