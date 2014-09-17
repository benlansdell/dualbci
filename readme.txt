==============================
BCI Target pursuit task -- GLM
==============================

MATLAB code for studying neural recording data from monkey as performs manual
and brain control target pursuit tasks. Data from Moritz lab. Orignally studied
by Charlie Matlack (cmatlack@uw.edu). 

Dependencies: makes use of functions in the stats toolbox.

In preprocess:
-Functions for importing data from .nev and .nsx files. 
-Functions for importing trial data from Labview (.mat) files. 
-Functions for smoothing torque data, and converting spike times to binned spikes
-generate_glm_data.m for simulating a GLM given input data and filters

In models:
-If the GLM has the form: E(y) = g(X\beta) where X is a data matrix and \beta is
 a vector of filter coefficients, this directory contains functions for taking raw 
 data and preparing data matrix to be input into a GLM, along with output vector y.

In eval:
-Functions to fit the filter coefficients beta
-Functions to predict spike trains given stimulus and a fit GLM
-Functions to plot filters of GLM
-Other plotting functions

In other: 
-Other people's GLM code

In old:
-Old code

In worksheets:
-Some scripts that make use of all this code

==========
How to use
==========

If MATLAB is started in the ./matlab directory then startup.m will
automatically add the above directories to the path. It will try to add
the Chronux functions by adding ~/matlab/chronux to the path, and will try to
add the .nev and .ns3 files to the path by looking in ./matlab/blackrock, but these
are not necessary.

See ./matlab/accessing_data.txt for information on the format of matlab,
labview and BlackRock files

===============
Version history
===============

0.1 -- contains a whole bunch of functions to compute correlation (corr_*.m) bw
firing rate and torque (cursor) position/vel/accel. Not needed, as superseded
by correlation_nev.m. Will be dropped from future versions.

