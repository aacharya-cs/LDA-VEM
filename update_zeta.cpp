/*********************************************************************
 * update_phi_cpp.cpp
 * @ Ayan Acharya, Date: Mar 31, 2012
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

void update_zeta_i(mxArray *retptr, double *sstopicwordptr, double *zetai, double *windexi, double *wcounti, const mxArray *model, const mxArray *data, const double *Esticks, const double *Elogbeta, const double *Elogmu, const double *smallphi, const int ndistWords, const int i, double *count2)
{
    int index, nK, t, k, k1, T, K1, K2, V, N, tempind1, tempind2, j, y, C2, Y, phase;
    double logsum, minval, val1, val2, val3, value, epsilon;
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
    nK      = (T+K2);
    epsilon = mxGetScalar(mxGetField(model,0,"epsilon"));
    
    if(phase==1)
    {// use the dual variable only in training phase; no dual variable in test phase
        dmu     = (double*)mxGetPr(mxGetField(model,0,"dmu"));
    }
    
    tmp    = mxCreateDoubleMatrix(ndistWords,nK,mxREAL);
    tmpptr = mxGetPr(tmp);
    
    for (j=0; j<ndistWords; j++) // loop over (distinct) words
    {
        logsum = 0;
        tmp1 = Malloc(double,nK);
        
        for (k=0; k<nK; k++)    // loop over third dimension
        {
            tempind1 = k + ((int)(windexi[j])-1)*(K1+K2);
            val1     = 0;
            
            if(phase==1 && (int)classlabels[i]>=1)   // use the dual variable only in training phase only when label is present; no dual variable in test phase;
            {
                for (y=0; y<Y; y++)
                {
                    val1 = val1 + dmu[y*N+i]*(r[k*Y+ (int)classlabels[i]-1] - r[k*Y+y]);
                }
                val1 =  val1/nwordspdoc[i];
            }
            
            if(k<T) // for unsupersvised topics
            {
                val2 = 0;
                for(k1=0; k1<K1; k1++)
                {
                    tempind2 = k1 + ((int)(windexi[j])-1)*(K1+K2);
                    val2     = val2 + smallphi[i+k*N+k1*N*T]*Elogbeta[tempind2];
                }
                
                /*if((int)count2[0]<3)
                    val3 = 0;
                else*/
                    val3 = Esticks[k*N+i];
                *(tmp1+k) = (val1 + val2 + val3) + log(1-epsilon);
            }
            else    // for supersvised topics
            {
                *(tmp1+k) = Elogmu[k*N+i] + Elogbeta[tempind1] + val1 + log(epsilon);
            }
            //if condition says when to ignore phi's
            if(phase==1) // only in training phase
            {
                if (k>=T && *(annotations+(k-T)*N+i)==0);
                //mexPrintf("%d %d %d ahh here\n", i, j, k);  // DSLDA
                else
                    logsum = log_sum(*(tmp1+k),logsum);
            }
            else
                logsum = log_sum(*(tmp1+k),logsum);
            
        }
        
        //mexPrintf("%d %d %f\n", i, j, logsum);
        
        // conversion from log space to real number
        for (k=0; k<nK; k++)
        {
            if(logsum - *(tmp1+k)>1000)
                tmpptr[k*ndistWords+j] = minval;
            if(logsum - *(tmp1+k)<1000)
                tmpptr[k*ndistWords+j] = exp(*(tmp1+k) - logsum) + minval;
            
            if(phase==1) // only in training phase
            {
                if (k>=T && *(annotations+(k-T)*N+i)==0)
                {
                    tmpptr[k*ndistWords+j] = 0;  // DSLDA
                    //mexPrintf("%d %d %d hey here!\n", i, j, k);
                }
                else; // do nothing
            }
        }
        free(tmp1);
    }
    
    // update sufficient statistics -- for both unsupervised and supervised topics
    for (k = 0; k < (K1+K2); k++)
    {
        for (j = 0; j < ndistWords; j++)
        {
            value = 0;
            if(k<K1)  // unsupervised topics
            {
                for (t = 0; t < T; t++) // loop over topics
                {
                    //mexPrintf("hey here1! %d %d %d %d %f\n", k, i, ndistWords, t, windexi[j]);
                    value += smallphi[i+t*N+k*N*T]*tmpptr[j+t*ndistWords]; // smallphi{itk1}*zeta{ijt}
                }
                value = wcounti[j]*value;
            }
            else   // supervised topics
            {
                value = wcounti[j]*tmpptr[j+(k-K1+T)*ndistWords];
            }
            index = k + ((int)windexi[j]-1)*(K1+K2);
            if(value<0)
                mexPrintf("error -- %d %d %f\n", k, j, tmpptr[j+(k-K1+T)*ndistWords]);
            sstopicwordptr[index] += value;
        }
    }
    mxSetCell (retptr, i, tmp);
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *model, *data;
    mxArray *zeta, *windex, *wcount;
    double  *zetai, *windexi, *wcounti, *Esticks, *Elogbeta, *Elogmu, *smallphi, *sstopicwordptr, *count2;
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
    Elogbeta     = (double*)mxGetPr(prhs[3]); // psi(\lambda)
    Elogmu       = (double*)mxGetPr(prhs[4]); // sum of psi(\lambda) over vocabulary
    count2       = (double*)mxGetPr(prhs[5]);
    
    plhs[0] = mxCreateCellMatrix(1,N);
    plhs[1] = mxCreateDoubleMatrix((K1+K2),V,mxREAL);
    sstopicwordptr = mxGetPr(plhs[1]);
    
    for (int i = 0; i < N; i++)
    {
        zetai   = (double*)mxGetPr(mxDuplicateArray (mxGetCell (zeta, i)));     // pointer to zeta{i}
        windexi = (double*)mxGetPr(mxDuplicateArray (mxGetCell (windex, i)));   // pointer to windex{i}
        wcounti = (double*)mxGetPr(mxDuplicateArray (mxGetCell (wcount, i)));   // pointer to wcount{i}
        ndistWords = mxGetM(mxDuplicateArray (mxGetCell (windex, i)));
        if(mxGetN(mxDuplicateArray (mxGetCell (windex, i)))>ndistWords)
            ndistWords = mxGetN(mxDuplicateArray (mxGetCell (windex, i)));
        //mexPrintf("nwords: %d\n", ndistWords);
        update_zeta_i(plhs[0], sstopicwordptr, zetai, windexi, wcounti, model, data, Esticks, Elogbeta, Elogmu, smallphi, ndistWords, i, count2);
    }
    
    return;
}
