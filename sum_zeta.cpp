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

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *tmp1, *tmp3, *tmp5, *tmp6;
    double *tmp2, *tmp4, *tmp7;
    int N; // number of instances
    int dimx, dimy;
    
    tmp1 = mxGetField(prhs[0],0,"zeta");
    N = mxGetNumberOfElements (tmp1);
    int K2 = (int)mxGetScalar(mxGetField(prhs[0],0,"K2"));
    int T = (int)mxGetScalar(mxGetField(prhs[0],0,"T"));
    int nK = (T+K2);
    
    plhs[0] = mxCreateNumericMatrix(N, nK, mxDOUBLE_CLASS, 0);
    tmp2 = mxGetPr(plhs[0]);

    tmp5 = mxGetField(prhs[1],0,"wcount");   
    
    for (int i = 0; i < N; i++)
    {
        tmp3 = mxDuplicateArray (mxGetCell (tmp1, i));
        tmp4 = mxGetPr(tmp3);
        dimx = (int)mxGetM(tmp3);
        dimy = (int)mxGetN(tmp3);        

        tmp6 = mxDuplicateArray (mxGetCell (tmp5, i));
        tmp7 = mxGetPr(tmp6);    

        //matlab's array index is column wise, so instead of accessing (j,k), access (k,j) when j indexes row and k indexes column       
        
        for (int k = 0; k < dimy; k++)
        {
            double value = 0;
            for(int j=0; j<dimx; j++)
            {
                value+= tmp7[j]*tmp4[k*dimx+j]; //w_{nm'}*zeta_{nm'k}
            }
            tmp2[k*N+i]+= value;
        }
    }
    return;
}


