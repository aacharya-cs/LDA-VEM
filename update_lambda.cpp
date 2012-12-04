/*********************************************************************
 * update_lambda.cpp
 * @ Ayan Acharya, Date: 11/8/2012
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
    const mxArray *model, *data;
    mxArray *zeta, *wcount, *windex;
    double value, *smallphi, *eta, *zetai, *lambda, *windexi, *wcounti;
    int i, j, k, t, v, nK, K1, K2, T, N, V, ndistWords;
    
    K1  = (int)mxGetScalar(mxGetField(prhs[0],0,"K1")); // number of latent topics
    K2  = (int)mxGetScalar(mxGetField(prhs[0],0,"K2")); // number of superised topics
    N   = (int)mxGetScalar(mxGetField(prhs[0],0,"N"));  // number of instances
    V   = (int)mxGetScalar(mxGetField(prhs[0],0,"V"));  // size   of vocabulary
    T   = (int)mxGetScalar(mxGetField(prhs[0],0,"T"));  // number of topics
    nK  = (K1+K2);
    
    model = prhs[0];
    data  = prhs[1];
    
    zeta       = mxGetField(model,0,"zeta");
    wcount     = mxGetField(data,0,"wcount");
    windex     = mxGetField(data,0,"windex");
    smallphi   = (double*)mxGetPr(mxGetField(model,0,"smallphi"));
    eta        = (double*)mxGetPr(mxGetField(model,0,"eta"));
    
    plhs[0] = mxCreateNumericMatrix(nK, V, mxDOUBLE_CLASS, 0); // lambda _{kv}
    lambda  = (double*)mxGetPr(plhs[0]);
    
    for (k = 0; k < nK; k++)        // loop over all topics
    {
        for (v = 0; v < V; v++)     // loop over vocabulary and intialize to eta[v]
            lambda[v*nK+k] = eta[v];
        
        for (i = 0; i < N; i++)     // loop over documents
        {
            zetai   = (double*)mxGetPr(mxDuplicateArray (mxGetCell (zeta, i)));     // zeta{i}..double array
            wcounti = (double*)mxGetPr(mxDuplicateArray (mxGetCell (wcount, i)));   // wcount{i}..integer array
            windexi = (double*)mxGetPr(mxDuplicateArray (mxGetCell (windex, i)));   // windex{i}..integer array
            ndistWords = mxGetM(mxDuplicateArray (mxGetCell (windex, i)));
            if(mxGetN(mxDuplicateArray (mxGetCell (windex, i)))>ndistWords)
                ndistWords = mxGetN(mxDuplicateArray (mxGetCell (windex, i)));
            
            for(j=0; j<ndistWords; j++)     // loop over (distinct) words
            {
                if(k<K1)                    // for latent topics
                {
                    value = 0;
                    for (t = 0; t < T; t++) // loop over topics
                    {
                        //mexPrintf("hey here1! %d %d %d %d %f\n", k, i, ndistWords, t, windexi[j]);
                        value += smallphi[i+t*N+k*N*T]*zetai[j+t*ndistWords]; // smallphi{itk1}*zeta{ijt}
                    }
                    lambda[(int)(windexi[j]-1)*nK+k] += wcounti[j]*value;
                }
                else                        // for supervised topics
                {                           
                    //mexPrintf("!!!!!error: word index zero!!!!\n");
                    value = zetai[j+(k-K1+T)*ndistWords];
                    lambda[(int)(windexi[j]-1)*nK+k] += wcounti[j]*value;
                }
            }
            //mexPrintf("hey here2! %d %d %d\n", k, i, ndistWords);
        }
    }
    return;
}