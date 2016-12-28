# Analysis of dual control BCI target pursuit task

Ben Lansdell, Ivana Milovanovic. 2016

MATLAB code for studying multi-electrode array recording data from monkey as performs manual, brain control and dual control target pursuit task. Data from [Moritz lab](http://depts.washington.edu/moritlab/), experiment setup by Charlie Matlack. 

## Dependencies:

A full working environment requires the following:

* `glmfit`, which is in the statistics toolbox.
* Mark Schmit's L1 toolbox ([here](http://www.cs.ubc.ca/~schmidtm/Software/thesis.html))
* [MVGC toolbox](http://users.sussex.ac.uk/~lionelb/MVGC/) 
* [GPFA toolbox](https://users.ece.cmu.edu/~byronyu/software.shtml)
* [Transfer entropy toolbox](https://code.google.com/archive/p/transfer-entropy-toolbox/) 
* [MySQL java driver](http://dev.mysql.com/downloads/connector/j/)

Depending on your usage, only a subset of these may be required. See usage section below.

## Code

### In preprocess:

preprocess\_\* functions generate a structure for use by functions in ./models

* Functions for importing data from .nev and .nsx files. 
* Functions for importing trial data from Labview (.mat) files. 
* Functions for smoothing torque data, and converting spike times to binned spikes
* generate_glm_data.m for simulating a GLM given input data and filters

### In models:

Functions to take a preprocess structure output from ./preprocess, and generate a structure with data in appropriate format for fitting a GLM/linear model/etc to.

If the GLM has the form: $E(y) = g(X\beta)$ where $X$ is a data matrix and $\beta$ is a vector of filter coefficients, this directory contains functions for taking raw data and preparing data matrix to be input into a GLM, along with output vector y.

### In fitting: 

Functions to perform MLE fit of GLM model

### In eval:

Functions to interpret results of MLE fitting

* Functions to predict spike trains given stimulus and a fit GLM
* Functions to plot filters of GLM
* Granger causality connectivity
* Other plotting functions

### In scripts:
* Scripts that generate plots used in paper

### In functions:

Various support functions

* Functions to import blackrock files
* Functions to compute correlations, etc
* Functions to save plots as .eps
* Various other things

## In sql:

The size of the dataset and number of analyses performed necessitated managing fits and results in a SQL database. Here are functions to add fits and other analyses to a SQL database. Please contact us directly if you'd like to access our database containing our stored analyses and results, or would like to setup something similar yourself. 

## Ok enough details, how do I use this:

The simplest use case is:

* Making plots presented in paper:
  See scripts in ./scripts

Beyond that, you'll need to set some paths:

1. If MATLAB is started in this directory then `startup.m` will automatically add the above directories to the path. If not, make sure `startup.m` is run to setup paths and toolboxes, etc. You'll need to open `startup.m` to check the paths to different toolboxes are set correctly.

2. The next simplest use case is to use the GLM code to fit/interpret different models:
```const = 'on';
nK_sp = 100; 
nK_pos = 100;
%Load test preprocessed data
pre = load('./testdata/test_preprocess_spline_short.mat');
data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
model = MLE_glmfit(data, const);
```

Beyond that, to further setup your environment to be able to redo all analyses performed in our paper:

* This requires acquiring the data. Once acquired its path can be added to `startup.m`. By default it looks for BlackRock .nev and .ns3 files in ./blackrock and ./labview. Please contact the authors for more information.

See ./accessing_data.txt for information on the format of matlab, labview and BlackRock files.
