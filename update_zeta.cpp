/*********************************************************************
 * update_zeta.cpp
 * @ Ayan Acharya, Date: Dec 6, 2012
 * input:
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

void update_zeta_i(mxArray *retptr, double *sstopicwordptr, double *ssfeatures, double *zetai, double *windexi, double *wcounti, const mxArray *model, const mxArray *data, const double *Esticks, const double *Eloglambda, const double *Eloggamma, const double *smallphi, const int ndistWords, const int i, double *count2, double *option)
{
    int index, nK1, nK2, t, k, k1, T, K1, K2, V, N, tempind1, tempind2, j, y, C2, Y, phase;
    double logsum, logsum1, logsum2, minval, val1, val2, val3, value, epsilon, valtemp;
    double *tmpptr, *tmp1, *log_beta, *dmu, *r, *annotations, *classlabels, *nwordspdoc;
    mxArray *tmp;
    
    minval      = mxGetScalar(mxGetField(model,0,"MINVALUE"));
    phase       = mxGetScalar(mxGetField(model,0,"phase"));
    classlabels = (double*)mxGetPr(mxGetField(data,0,"classlabels"));
    nwordspdoc  = (double*)mxGetPr(mxGetField(data,0,"nwordspdoc"));
    annotations = (double*)mxGetPr(mxGetField(data,0,"annotations"));
    
    r       = mxGetPr(mxGetField(model,0,"r"));
    V       = mxGetScalar(mxGetField(model,0,"V"));
    N       = mxGetScalar(mxGetField(model,0,"N"));
    C2      = mxGetScalar(mxGetField(model,0,"C2"));
    Y       = mxGetScalar(mxGetField(data,0,"Y"));
    T       = mxGetScalar(mxGetField(model,0,"T"));
    K1      = mxGetScalar(mxGetField(model,0,"K1"));
    K2      = mxGetScalar(mxGetField(model,0,"K2"));
    epsilon = mxGetScalar(mxGetField(model,0,"epsilon"));
    
    if((int)option[0]==1)
    {
        nK1 = (T+K2);
        nK2 = (K1+K2);
    }
    else
    {
        nK1 = T;
        nK2 = K1;
    }
    
    if(phase==1)
    {
        // use the dual variable only in training phase; no dual variable in test phase
        dmu     = (double*)mxGetPr(mxGetField(model,0,"dmu"));
    }
    
    tmp    = mxCreateDoubleMatrix(ndistWords,nK1,mxREAL);
    tmpptr = mxGetPr(tmp);
    
    for (j=0; j<ndistWords; j++)  // loop over (distinct) words
    {
        logsum  = 0;
        logsum1 = 0;
        logsum2 = 0;
        
        tmp1 = Malloc(double,nK1);
        
        for (k=0; k<nK1; k++)     // loop over third dimension
        {
            tempind1 = k + ((int)(windexi[j])-1)*nK2;
            val1     = 0;
            
            ////////////////////////////////////////////////////////////////////////////////////////////
            // terms from document level supervision
            if(phase==1 && (int)classlabels[i]>=1)   // use the dual variable only in training phase only when label is present; no dual variable in test phase;
            {
                if(k<T)          // for unsupervised topics
                {
                    for (y=0; y<Y; y++)
                    {
                        valtemp = 0;
                        for (k1=0; k1<K1; k1++)
                        {
                            valtemp = valtemp + (r[k1*Y+ (int)classlabels[i]-1] - r[k1*Y+y])*smallphi[i+k*N+k1*N*T];
                        }
                        val1 = val1 + dmu[i+y*N]*valtemp;
                    }
                }
                else   // for NPDSLDA
                    if((int)option[0]==1)// for supervised topics
                    {
                    for (y=0; y<Y; y++)
                    {
                        val1 = val1 + dmu[i+y*N]*(r[(int)classlabels[i]-1 + (k-T+K1)*Y] - r[y + (k-T+K1)*Y]);
                    }
                    val1 =  val1/nwordspdoc[i];
                    }
            }
            
            ////////////////////////////////////////////////////////////////////////////////////////////
            // for other terms
            if(k<T) // for unsupersvised topics
            {
                val2 = 0;
                for(k1=0; k1<K1; k1++)
                {
                    tempind2 = k1 + ((int)(windexi[j])-1)*nK2;
                    val2     = val2 + smallphi[i+k*N+k1*N*T]*Eloglambda[tempind2];
                }
                val3 = Esticks[i+k*N];
                *(tmp1+k) = val1 + val2 + val3;
            }
            else    // for supersvised topics
                if((int)option[0]==1)  // for NPDSLDA
                {
                *(tmp1+k) = Eloggamma[i+k*N] + Eloglambda[tempind1] + val1;
                }
            
            
            if(phase==1) // only in training phase; no need to have any clause for test phase
            {
                if((int)option[0]==1) // NPDSLDA
                {
                    if (k<T)    // unsupervised topics
                        logsum1 = log_sum(*(tmp1+k),logsum1);
                    else        // supervised topics
                    {
                        if(*(annotations+(k-T)*N+i)==0);  //if condition says when to ignore phi's
                        else
                            logsum2 = log_sum(*(tmp1+k),logsum2);
                    }
                }
                else                  // NPLDA
                    logsum = log_sum(*(tmp1+k),logsum);
            }
        }
        
        // conversion from log space to real number
        for (k=0; k<nK1; k++)
        {
            if((int)option[0]==1) // NPDSLDA
            {
                if(k<T)     // unsupervised topics
                {
                    if(logsum1 - *(tmp1+k)>10000)
                        tmpptr[k*ndistWords+j] = minval;
                    if(logsum1 - *(tmp1+k)<10000)
                        tmpptr[k*ndistWords+j] = (1-epsilon)*exp(*(tmp1+k) - logsum1) + minval;
                }
                else        // supervised topics
                {
                    if(logsum2 - *(tmp1+k)>10000)
                        tmpptr[k*ndistWords+j] = minval;
                    if(logsum2 - *(tmp1+k)<10000)
                        tmpptr[k*ndistWords+j] = epsilon*exp(*(tmp1+k) - logsum2) + minval;
                    if (phase==1 && k>=T && *(annotations+(k-T)*N+i)==0) // only in training phase
                    {
                        tmpptr[k*ndistWords+j] = 0;  // DSLDA
                        //mexPrintf("%d %d %d hey here!\n", i, j, k);
                    }
                }
            }
            else  // NPLDA
            {
                if(logsum - *(tmp1+k)>10000)
                    tmpptr[k*ndistWords+j] = minval;
                if(logsum - *(tmp1+k)<10000)
                    tmpptr[k*ndistWords+j] = exp(*(tmp1+k) - logsum) + minval;
            }
        }
        
        free(tmp1);
    }
    
    // update sufficient statistics -- for both unsupervised and supervised topics
    for (k = 0; k < nK2; k++)
    {
        for (j = 0; j < ndistWords; j++)
        {
            value = 0;
            if(k<K1)  // unsupervised topics
            {
                for (t = 0; t < T; t++) // loop over topics
                {
                    //mexPrintf("hey here1! %d %d %d %d %f\n", k, i, ndistWords, t, windexi[j]);
                    value += smallphi[i+t*N+k*N*T]*tmpptr[j+t*ndistWords]; // smallphi{ntk1}*zeta{nmt}
                }
                value = wcounti[j]*value;
            }
            else    // for NPDSLDA
                if((int)option[0]==1)// supervised topics
                {
                value = wcounti[j]*tmpptr[j+(k-K1+T)*ndistWords]; // zeta{nmk2}
                }
            index = k + ((int)windexi[j]-1)*nK2;
            sstopicwordptr[index] += value;
            ssfeatures[i+k*N] += value;
        }
    }
    
    mxSetCell (retptr, i, tmp);
    
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *model, *data;
    mxArray *zeta, *windex, *wcount;
    double  *zetai, *windexi, *wcounti, *Esticks, *Eloglambda, *Eloggamma;
    double  *smallphi, *sstopicwordptr, *ssfeatures, *count2, *option;
    int N, ndistWords, K1, K2, V;
    
    model = prhs[0];
    data  = prhs[1];
    N     = (int)mxGetScalar(mxGetField(model,0,"N"));
    K1    = (int)mxGetScalar(mxGetField(model,0,"K1"));
    K2    = (int)mxGetScalar(mxGetField(model,0,"K2"));
    V     = (int)mxGetScalar(mxGetField(model,0,"V"));
    
    zeta         = mxGetField(model,0,"zeta");
    windex       = mxGetField(data,0,"windex");
    wcount       = mxGetField(data,0,"wcount");
    smallphi     = (double*)mxGetPr(mxGetField(model,0,"smallphi"));
    Esticks      = (double*)mxGetPr(prhs[2]); // E_{log(\pi)}
    Eloglambda   = (double*)mxGetPr(prhs[3]); // psi(\lambda)
    Eloggamma    = (double*)mxGetPr(prhs[4]); // sum of psi(\lambda) over vocabulary
    count2       = (double*)mxGetPr(prhs[5]); // counter for controlling initialization
    option       = (double*)mxGetPr(prhs[6]); // option=1, NPDSLDA; (int)option[0]==0, NPLDA
    
    plhs[0] = mxCreateCellMatrix(1,N);
    if((int)option[0]==1)
    {
        plhs[1] = mxCreateDoubleMatrix((K1+K2),V,mxREAL);
        plhs[2] = mxCreateDoubleMatrix(N,(K1+K2),mxREAL);
    }
    else
    {
        plhs[1] = mxCreateDoubleMatrix(K1,V,mxREAL);
        plhs[2] = mxCreateDoubleMatrix(N,K1,mxREAL);
    }
    sstopicwordptr = mxGetPr(plhs[1]);
    ssfeatures = mxGetPr(plhs[2]);
    
    for (int i = 0; i < N; i++)
    {
        zetai   = (double*)mxGetPr(mxDuplicateArray (mxGetCell (zeta, i)));     // pointer to zeta{i}
        windexi = (double*)mxGetPr(mxDuplicateArray (mxGetCell (windex, i)));   // pointer to windex{i}
        wcounti = (double*)mxGetPr(mxDuplicateArray (mxGetCell (wcount, i)));   // pointer to wcount{i}
        ndistWords = mxGetM(mxDuplicateArray (mxGetCell (windex, i)));
        if(mxGetN(mxDuplicateArray (mxGetCell (windex, i)))>ndistWords)
            ndistWords = mxGetN(mxDuplicateArray (mxGetCell (windex, i)));
        //mexPrintf("hey here! %d\n", i);
        update_zeta_i(plhs[0], sstopicwordptr, ssfeatures, zetai, windexi, wcounti, model, data, Esticks, Eloglambda, Eloggamma, smallphi, ndistWords, i, count2, option);
    }
    return;
}
