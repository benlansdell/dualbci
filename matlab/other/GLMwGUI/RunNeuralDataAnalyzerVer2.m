clear all
close all
clc

% adding paths
basedir = pwd;
addpath([basedir '/tools_mex/']);
addpath([basedir '/tools_revcorrelations/']);
addpath([basedir '/tools_newglm/']);
addpath([basedir '/tools_subfig/']);

% compiling mex files
cd tools/mex

% we compile different way depending on the operating system.
if exist('crosscorr_mex.mexmaci64', 'file') == 0    
    % some works with only c++ or some works only with c
    mex crosscorr_mex.cpp
    mex crosscorr_mex.c
end

cd ..

NeuralDataAnalyzer_v2