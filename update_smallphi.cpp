/*********************************************************************
 * update_smallphi.cpp
 * @ Ayan Acharya, Date: 11/8/2012
 * input: 1: model, 2: data, 3: psi(1,gamma), 4: option, 5: phase, 6: annotation,
 ********************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <matrix.h>
#include <mex.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

double log_sum(double log_a, double log_b)
{
    double v;
    if (log_a == 0) return(log_b);
    if (log_b == 0) return(log_a);
    if (log_a < log_b)
        v = log_b+log(1 + exp(log_a-log_b));
    else
        v = log_a+log(1 + exp(log_b-log_a));
    return(v);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *model, *data;
    mxArray *windex, *wcount, *zeta;
    double  *Esticks, *Elogbeta, *dmu, *r, *classlabels, *sumzeta, *nwordspdoc;
    double  *smallphi, *zetai, *tmp, *windexi, *wcounti, *count2;
    double  minval, value1, value2, value3;
    int     *dims, N, K1, K2, T, V, i, t, k1, w, y, Y, tempind, ndistWords, phase;
    
    model = prhs[0];
    data  = prhs[1];
    
    T      = (int)mxGetScalar(mxGetField(model,0,"T"));
    K1     = (int)mxGetScalar(mxGetField(model,0,"K1"));
    K2     = (int)mxGetScalar(mxGetField(model,0,"K2"));
    N      = (int)mxGetScalar(mxGetField(model,0,"N"));
    V      = (int)mxGetScalar(mxGetField(model,0,"V"));
    phase  = (int)mxGetScalar(mxGetField(model,0,"phase"));
    minval = (double)mxGetScalar(mxGetField(model,0,"MINVALUE"));
    
    windex        = mxGetField(data,0,"windex");
    wcount        = mxGetField(data,0,"wcount");
    zeta          = mxGetField(model,0,"zeta");
    sumzeta       = (double*)mxGetPr(mxGetField(model,0,"sumzeta"));
    nwordspdoc    = (double*)mxGetPr(mxGetField(data,0,"nwordspdoc"));
    classlabels   = (double*)mxGetPr(mxGetField(data,0,"classlabels"));
    Esticks       = (double*)mxGetPr(prhs[2]); // Eq[log(\beta)]
    Elogbeta      = (double*)mxGetPr(prhs[3]); // psi(\lambda) - psi(\sum_{v=1}^{V} \lambda)
    count2        = (double*)mxGetPr(prhs[4]); // counter for controlling initialization
    
    if(phase==1)
    {// use the dual variable only in training phase; no dual variable in test phase
        dmu     = (double*)mxGetPr(mxGetField(model,0,"dmu"));
    }
    
    dims      = Malloc(int,3);
    *dims     = N;
    *(dims+1) = T;
    *(dims+2) = K1;
    plhs[0]   = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    smallphi  = (double*)mxGetPr(plhs[0]); // smallphi_{itk1}
    free(dims);
    
    for (i=0; i<N; i++)                    // loop over documents
    {
        zetai   = (double*)mxGetPr(mxDuplicateArray (mxGetCell (zeta, i)));     // zeta_{i} double array
        windexi = (double*)mxGetPr(mxDuplicateArray (mxGetCell (windex, i)));   // pointer to windex{i}
        wcounti = (double*)mxGetPr(mxDuplicateArray (mxGetCell (wcount, i)));   // pointer to wcount{i}
        
        ndistWords = mxGetM(mxDuplicateArray (mxGetCell (windex, i)));
        if(mxGetN(mxDuplicateArray (mxGetCell (windex, i)))>ndistWords)
            ndistWords = mxGetN(mxDuplicateArray (mxGetCell (windex, i)));
        //mexPrintf("ok till here1 %d %d\n", i, ndistWords);
        
        for (t=0; t<T; t++)           // loop over topics (lower level truncation)
        {
            double logsum = 0;
            tmp = Malloc(double,K1);
            
            for (k1=0; k1<K1; k1++)   // loop over topics (higher level truncation)
            {
                value1 = 0;
                value3 = 0;
                if(phase==1 && (int)classlabels[i]>=1)   // use the dual variable only in training phase only when label is present; no dual variable in test phase;
                {
                    for (y=0; y<Y; y++)
                    {
                        value3 = value3 + dmu[y*N+i]*(r[k1*Y+ (int)classlabels[i]-1] - r[k1*Y+y]);
                    }
                    value3 =  value3*sumzeta[i+t*N]/nwordspdoc[i];
                }
                for (w=0; w<ndistWords; w++) // loop over (distinct) words
                {
                    tempind = k1 + ((int)(windexi[w])-1)*K1;
                    value1  = value1 + wcounti[w]*zetai[w+t*ndistWords]*Elogbeta[tempind];
                }
                /*if((int)count2[0]<3)
                 * value2  = 0;
                 * else*/
                value2  = Esticks[k1];
                *(tmp+k1) = (value1 + value2 + value3);
                logsum    = log_sum(*(tmp+k1),logsum);
            }
            // conversion from log space to real number
            for (k1=0; k1<K1; k1++)
            {
                if(logsum - *(tmp+k1)>1000)
                    smallphi[i+t*N+k1*N*T] = minval;
                if(logsum - *(tmp+k1)<1000)
                    smallphi[i+t*N+k1*N*T] = exp(*(tmp+k1)-logsum)+minval;
            }
            free(tmp);
        }
    }
    
    return;
}
