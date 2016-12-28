# Analysis of dual control BCI target pursuit task

Ben Lansdell, Ivana Milovanovic. 2016

MATLAB code for studying multi-electrode array recording data from monkey as performs manual, brain control and dual control target pursuit task. Data from [Moritz lab](http://depts.washington.edu/moritlab/), experiment setup by Charlie Matlack. 

### Dependencies:
* Makes use of functions in the stats toolbox.

### In preprocess:
* Functions for importing data from .nev and .nsx files. 
* Functions for importing trial data from Labview (.mat) files. 
* Functions for smoothing torque data, and converting spike times to binned spikes
* generate_glm_data.m for simulating a GLM given input data and filters

### In models:
* If the GLM has the form: E(y) = g(X\beta) where X is a data matrix and \beta is
 a vector of filter coefficients, this directory contains functions for taking raw 
 data and preparing data matrix to be input into a GLM, along with output vector y.

### In eval:
* Functions to fit the filter coefficients beta
* Functions to predict spike trains given stimulus and a fit GLM
* Functions to plot filters of GLM
* Other plotting functions

### In functions:
* Functions to import blackrock files
* Functions to compute correlations, etc
* Functions to save plots as .eps
* Some other things

### In scripts:
* Scripts that generate plots used in paper

## How to use

If MATLAB is started in this directory then startup.m will automatically add the above directories to the path. It will try to add the Chronux functions by adding ~/matlab/chronux to the path, and will try to add the .nev and .ns3 files to the path by looking in ./matlab/blackrock, but these are not necessary.

See ./accessing_data.txt for information on the format of matlab, labview and BlackRock files
