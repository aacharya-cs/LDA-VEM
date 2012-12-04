/*********************************************************************
 * update_beta_cpp.cpp
 * @ Ayan Acharya, Date: Mar 31, 2012
 ********************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <matrix.h>
#include <mex.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *tmp1, *tmp2, *tmp3;
    int     dimx, dimy, ii, jj, i, j, nK;
    
    //mexPrintf("ok till here 1\n");
    
    int N  = (int)mxGetScalar(mxGetField(prhs[0],0,"N"));
    int V  = (int)mxGetScalar(mxGetField(prhs[0],0,"V"));
    int option  = (int)mxGetScalar(mxGetField(prhs[0],0,"option"));
    
    nK = (int)mxGetScalar(mxGetField(prhs[0],0,"K"));
    
    //mexPrintf("ok till here 2\n");
    
    tmp1 = mxGetPr(mxGetField(prhs[0],0,"ss_topicword"));
    tmp2 = mxGetPr(mxGetField(prhs[0],0,"ss_topic"));
    plhs[0] = mxCreateNumericMatrix(nK, V, mxDOUBLE_CLASS, 0);
    tmp3    = mxGetPr(plhs[0]);
    
    for (ii = 0; ii < nK; ii++)             // loop over topics
    {
        //mexPrintf("topic: %d\n",ii);
        for (jj = 0; jj < V; jj++)          // loop over words in a vocabulary
        {
            if (tmp1[ii+jj*nK] > 0)
                tmp3[ii+jj*nK] = log(tmp1[ii+jj*nK]) - log(tmp2[ii]);
            else
                tmp3[ii+jj*nK] = -100;
        }
    }
    
    return;
}


