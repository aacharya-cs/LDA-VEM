/*********************************************************************
 * sum_phi.cpp
 * @ Ayan Acharya, Date: Mar 31, 2012
 ********************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *tmp1, *tmp3, *tmp5, *tmp6;
    double *tmp2, *tmp4, *tmp7;
    mwSize dimx, dimy, dimz, nK;
    
    tmp1 = prhs[0];
    dimx = mxGetDimensions(tmp1)[0];
    dimy = mxGetDimensions(tmp1)[1];
    dimz = mxGetDimensions(tmp1)[2];   
    tmp2 = (double*)mxGetPr(tmp1);
    
    //plhs[0] = mxCreateNumericMatrix(dimx, dimy, dimz, mxDOUBLE_CLASS, 0); // lambda _{kv}
    //tmp2 = mxGetPr(plhs[0]);
    
    for (int i = 0; i < dimx; i++)  // loop over number of latent topics
    {
        for (int j = 0; j < dimy; j++)  // loop over vocabulary
        {
            for (int k = 0; k < dimz; k++) // loop over documents
            {                
                mexPrintf("%d %d %d %f\n", i+1, j+1, k+1, tmp2[i+j*dimx+k*dimx*dimy]);
            }
        }
    }
    return;
}


