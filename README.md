# argo_fstmr_oxygen


This repository houses code for an revision of https://arxiv.org/abs/2211.04012, "A functional regression model for heterogeneous BioGeoChemical Argo data in the Southern Ocean"; Moritz Korte-Stapff, Drew Yarger, Stilian Stoev, and Tailen Hsing. A previous version of the code is here: https://github.com/dyarger/argo_clustering_cokriging.


There are two main components of our code. The first is the code that implements our methodology in a package, based on the folder "fstmr." Each file contains the relevant code for different parts of the methodology. There are three classes: one_var (for univariate data), mult_var_ind (for multivariate data with spatial independence), and mult_var_spat (our main methodology). The second part of the code is in paper/code, which contains the code for the data analysis, simulations, images, and processing results that are described in the paper. 

Major sources of data:

1. SOCCOM Data: https://library.ucsd.edu/dc/
object/bb9362707g

2. Argo Data: http://doi.org/10.17882/42182#85023

Minor sources of data:

1. RG climatology: https://sio-argo.ucsd.edu/RG_Climatology.html, ``2004-2018 1/6 x 1/6 degree mean field [t=0]'' available for download at https://sio-argo.ucsd.edu/gilson/argo_climatology/RG_ArgoClim_33pfit_2019_mean.nc.gz

2. SOSE grid: http://sose.ucsd.edu/BSOSE6_iter133_solution.html, with a direct download at http://sose.ucsd.edu/SO6/SETUP/grid.nc

3. Sea ice extent: National Snow and Ice Data Center, https://masie_web.apps.nsidc.org/pub/DATASETS/NOAA/G02135/south/monthly/geotiff/09_Sep/S_201409_concentration_v3.0.tif


