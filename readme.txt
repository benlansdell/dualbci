=======================
BCI Target pursuit task
=======================

MATLAB code for studying neural recording data from monkey as performs manual
and brain control target pursuit tasks. Data from Moritz lab. Orignally studied
by Charlie Matlack (cmatlack@uw.edu). 

Dependencies: Some of the LFP code assumes the package Chronux is available in
the path. (http://chronux.org)

In functions/nev:
- correlation_nev.m performs cross- and auto-correlation analysis to identity
  time lags and units that are linearly correlated with task
- regression_nev.m fits a linear filter between firing rate and torque
In functions/labview
- import_trials.m imports all trial information from a labview matlab file
- functions to plot trial info, etc
In functions/Tuning
- Charlie's functions to compute Tuning: cross-correlation and linear regression
In functions/LogProcess
- Charlie's functions to process labview files and save as matlab files,
  determine the 'mapping' (brain control, manual control, etc) used.
In functions/NPMK
- BlackRocks's Neural Processing MATLAB kit for reading .nev and .ns3 files into
  Matlab.

If MATLAB is started in the ./matlab directory then startup.m will
automatically add the above ./function directories. It will try to add
the Chronux functions by adding ~/matlab/chronux to the path, and will try to
add the .nev and .ns3 files to the path by looking in ./matlab/blackrock.

See ./matlab/accessing_data.txt for information on the format of matlab,
labview and BlackRock files

===============
Version history
===============

0.1 -- contains a whole bunch of functions to compute correlation (corr_*.m) bw
firing rate and torque (cursor) position/vel/accel. Not needed, as superseded
by correlation_nev.m. Will be dropped from future versions.

