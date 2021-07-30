This version is reliant upon estimation of technical white noise from timeseries collected from regions of the dish not containing cells, i.e. background fluorescence. 
The code uses a global Gaussian Process to infer the white noise variance from the entire dataset using joint log-likelihood, following which each
timeseries representing a cell is further analysed using a periodic and aperiodic covariance and a false discovery rate test is used to detect periodic timeseries.

To use this please cite 
The advantages of this approach for fluorescence imaging which contains auto-fluorescence are described in Manning et al (2019) "Quantitative single-cell live imaging links HES5 dynamics with cell-state and fate in murine neurogenesis
Nat Comms 10:2835. 
