This version is reliant upon estimation of technical white noise from timeseries collected from regions of the dish not containing cells, i.e. background fluorescence. 
The code uses a global Gaussian Process to infer the white noise variance from the entire dataset using joint log-likelihood, following which each
timeseries representing a cell is further analysed using a periodic and aperiodic covariance and a false discovery rate test is used to detect periodic timeseries.

To use this please cite the following publications:

Phillips et al. (2017) "Identifying stochastic oscillations in single-cell live imaging time series using Gaussian processes" PLoS Comput Biol 13(5): e1005479. Here you can find the the covariance models and the overall approach.

Manning et al. (2019) "Quantitative single-cell live imaging links HES5 dynamics with cell-state and fate in murine neurogenesis" Nat Comms 10:2835. Here you can find a discussion of the advantages of the joint likelihood approach for fluorescence imaging which contains auto-fluorescence.   

C Rasmussen and N Hannes (2010) "Gaussian processes for machine learning (GPML) toolbox" Journal of Machine Learning Research, 11:3011 3015.  Here you can find the information about the GPML toolbox. The copy shared with this repository contains additional routines. Please unarchive GPML.zip before running.



