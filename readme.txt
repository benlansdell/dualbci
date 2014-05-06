=========================
=BCI Target pursuit task=
=========================

MATLAB code for studying neural recording data from monkey as performs manual and brain control target pursuit tasks. Data from Moritz lab. Orignally studied by Charlie Matlack (cmatlack@uw.edu). 

Features:
- correlation_nev.m performs cross- and auto-correlation analysis to identity time lags and units that are linearly correlated with task
- regression_nev.m fits a linear filter between firing rate and torque
- import_trials.m imports all trial information from a labview matlab file

See ./matlab/accessing_data.txt for information on the format of matlab, labview and BlackRock files

If MATLAB is started in the ./matlab directory then startup.m will automatically add the ./function directories 

Version history:

0.1 -- contains a whole bunch of functions to compute correlation bw firing rate and torque(cursor) position/speed/accel. Not needed, to be dropped from future versions.
