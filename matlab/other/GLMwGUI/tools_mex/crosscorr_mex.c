/*
 * proj = crosscorr_mex(stim,feature)
 *
 * FeatSelect_mex.c
 * This is a MEX file to perform the cross-correlation of 2 vectors
 *
 * Author - SangWook Lee
 * Fairhall Lab.
 * University of Washington
 * Mar. 2014.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], /* Output variables */
        int nrhs, const mxArray *prhs[])     /* Input variables */
{
    
    double *temp, *proj, *feature, *stim;
    int i, ii, Nstim, Mstim, Nfeat, Mfeat, Mspt, Nspt;
    
    Mstim = mxGetM(prhs[0]); 
    
    Nstim = mxGetN(prhs[0]); /* Get a dimension of stim */
    
    Mfeat = mxGetM(prhs[1]);
    
    Nfeat = mxGetN(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(Mstim, Nstim, mxREAL);
    
    proj = mxGetPr(plhs[0]);
    
    feature = mxGetPr(prhs[1]);
    
    stim = mxGetPr(prhs[0]);
    
    for (i=0;i<(Nstim-Nfeat+1);i++){
        proj[i] = 0;
        for (ii = 0; ii < Nfeat; ii++){
                proj[i] = proj[i] + (feature[ii]*stim[i+ii]);
            }
    }
    
    for (i=(Nstim-Nfeat+1);i<Nstim;i++){
        proj[i] = 0;
        for (ii = 0; ii < (Nstim-i); ii++){
                proj[i] = proj[i] + (feature[ii]*stim[i+ii]);
            }
    }
}