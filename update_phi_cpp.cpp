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

void update_phi_n(mxArray *retptr, double *sstopicwordptr, double *sstopicptr, const mxArray *phi_n, const double *windexn, const double *wcountn, const mxArray *model, const mxArray *data, const double *psigammaptr, const int n, const int phase, const double *annotations, int option)
{
    int ndistWords, nK, k1, k2, V, N, tempind, i, j, y, C2, Y; // number of words, maximum number of topics, maximum number of observed topics
    double logsum, minval, val, epsilon;
    mxArray *tmp;
    double *tmpptr, *tmp1, *log_beta, *mu, *eta, *classlabels;
    double *nwordspdoc = mxGetPr(mxGetField(data,0,"nwordspdoc"));
    
    minval = mxGetScalar(mxGetField(model,0,"MINVALUE"));
    ndistWords = mxGetM(phi_n);
    nK     = mxGetN(phi_n);
    tmp    = mxCreateDoubleMatrix(ndistWords,nK,mxREAL);
    tmpptr = mxGetPr(tmp);
    log_beta = mxGetPr(mxGetField(model,0,"log_beta"));
    
    V      = mxGetScalar(mxGetField(model,0,"V"));
    N      = mxGetScalar(mxGetField(model,0,"N"));
    
    //mexPrintf("till here ok1\n");
    if(option>=3)
    {
        classlabels = mxGetPr(mxGetField(data,0,"classlabels"));
        if(phase==1)
        {// use the dual variable only in training phase; no dual variable in test phase
            mu     = mxGetPr(mxGetField(model,0,"mu"));
        }
        eta     = mxGetPr(mxGetField(model,0,"eta"));
        C2 = (int)mxGetScalar(mxGetField(model,0,"C2"));
        Y  = (int)mxGetScalar(mxGetField(model,0,"Y"));
        if(option>=4)
        {
            k1 = mxGetScalar(mxGetField(model,0,"k1"));
            k2 = mxGetScalar(mxGetField(model,0,"k2"));
            epsilon = mxGetScalar(mxGetField(model,0,"epsilon"));
            //mexPrintf("\nY: %d\n",Y);
        }
    }
    
    int lowlimit, uplimit;
    if(option==5 || option==8) // DSLDA-NSLT1, DSLDA-NSLT2
    {
        lowlimit = k1 + ((int)classlabels[n]-1)*(k2/Y);
        uplimit  = k1 + ((int)classlabels[n])*(k2/Y)-1;
        //mexPrintf("%d %d %d %d %d %d\n", n, k1, k2, lowlimit, uplimit, (int)classlabels[n]);
    }
    if(option==6) // DSLDA-NSLT-Ayan
    {
        lowlimit = k1 + ((int)classlabels[n]-1)*k2;
        uplimit  = k1 + ((int)classlabels[n])*k2-1;
        //mexPrintf("%d %d %d %d %d %d\n", n, k1, k2, lowlimit, uplimit, (int)classlabels[n]);
    }
    //mexPrintf("\t %d acajcjac %d",n, ndistWords);
    for (i=0; i<ndistWords; i++)
    {
        logsum = 0;
        tmp1 = Malloc(double,nK);
        //mexPrintf("%d %d %d %d till here ok2\n", n, i, windexn[i], ndistWords);
        
        for (j=0; j<nK; j++)
        {
            tempind = j + ((int)(windexn[i])-1)*nK;
            //mexPrintf("\n%d %d %d %d %f\n",n, i,j, nK, log_beta[tempind]);
            val = 0;
            if(phase==1 && (int)classlabels[n]>=1)   // use the dual variable only in training phase only when label is present; no dual variable in test phase;
            {
                for (y=0; y<Y; y++)
                    val = val + mu[y*N+n]*(eta[j*Y+ (int)classlabels[n]-1] - eta[j*Y+y]);
                val =  val/nwordspdoc[n];
            }
            
            //mexPrintf("\n%d %d %d %d %f\n",n, i,j, nK, log_beta[tempind]);
            if(option==3)
                *(tmp1+j) = psigammaptr[j*N+n] + log_beta[tempind] + val;  // access (j,i) th element from gamma
            
            if(option>=4)
            {
                if(j<k1)
                    *(tmp1+j) = psigammaptr[j*N+n] + log_beta[tempind] + log(epsilon) + val;  // access (j,i) th element from gamma
                else
                    *(tmp1+j) = psigammaptr[j*N+n] + log_beta[tempind] + log(1-epsilon) + val;  // access (j,i) th element from gamma
            }
            
            if(option==1 || option==2)
                *(tmp1+j) = psigammaptr[j*N+n] + log_beta[tempind];  // access (j,i) th element from gamma
            
            //if condition says when to ignore phi's
            if(phase==1 && option!=3) // only in training phase and not for medLDA
            {
                if(option==2 && *(annotations+j*N+n)==0);
                else if (option==4 && j<k1 && *(annotations+j*N+n)==0);  // DSLDA
                else if (option==5 && ((j<k1 && *(annotations+j*N+n)==0) || (j>=k1 && !(j>=lowlimit && j<=uplimit)))); // DSLDA-NSLT1
                else if (option==8 && ((j<k1)|| (j>=k1 && !(j>=lowlimit && j<=uplimit)))); // DSLDA-NSLT2
                else if (option==6 && ((j<k1 && *(annotations+j*N+n)==0) || (j>=k1 && !(j>=lowlimit && j<=uplimit)))); // DSLDA-NSLT-Ayan
                else if (option==7 && (j<k1 && *(annotations+j*N+n)==0)); // DSLDA-OSST
                else
                    logsum    = log_sum(*(tmp1+j),logsum);
                /*if(phase==1 && ((option==2 || (option==4 && j<k1)) && *(annotations+j*N+n)==0 || (option==5 && ((j<k1 && *(annotations+j*N+n)==0)||(!(j>=lowlimit && j<=uplimit)))) || )); // only in training phase*/
            }
            else
                logsum    = log_sum(*(tmp1+j),logsum);
        }
        
        // conversion from log space to real number
        for (int j=0; j<nK; j++)
        {
            if(logsum - *(tmp1+j)>50)
                tmpptr[j*ndistWords+i] = minval;                       //(j,i) th element
            if(logsum - *(tmp1+j)<50)
                tmpptr[j*ndistWords+i] = exp(*(tmp1+j)-logsum)+minval; //(j,i) th element
            
            if(phase==1 && option!=3) // only in training phase and not for medLDA or DSLDA-OSST
            {
                if(option==2 && *(annotations+j*N+n)==0)
                    tmpptr[j*ndistWords+i] = 0;
                else if (option==4 && j<k1 && *(annotations+j*N+n)==0)
                    tmpptr[j*ndistWords+i] = 0;  // DSLDA
                else if (option==7 && (j<k1 && *(annotations+j*N+n)==0))
                    tmpptr[j*ndistWords+i] = 0; // DSLDA-OSST; zero out in training phase; don't zero out latent topics, will handle that through epsilon
                else if (option==5 && ((j<k1 && *(annotations+j*N+n)==0) || (j>=k1 && !(j>=lowlimit && j<=uplimit))))
                {
                    //mexPrintf("%d\t %d\t %f\n", n, j, *(annotations+j*N+n));
                    tmpptr[j*ndistWords+i] = 0; // DSLDA-NSLT-Ray
                }
                else if (option==8 && ((j<k1) || (j>=k1 && !(j>=lowlimit && j<=uplimit)))) // DSLDA-NSLT2
                {
                    //mexPrintf("%d\t %d\t %f\n", n, j, *(annotations+j*N+n));
                    tmpptr[j*ndistWords+i] = 0; // DSLDA-NSLT-Ray
                }
                else if (option==6 && ((j<k1 && *(annotations+j*N+n)==0) || (j>=k1 && !(j>=lowlimit && j<=uplimit))))
                    tmpptr[j*ndistWords+i] = 0; // DSLDA-NSLT-Ayan
                else; // do nothing
                /*if(phase==1 && ((option==2 || (option==4 && j<k1)) && *(annotations+j*N+n)==0 || (option==5 && ((j<k1 && *(annotations+j*N+n)==0)||(!(j>=lowlimit && j<=uplimit))))))
                 * tmpptr[j*ndistWords+i] = 0; */
            }
            
            
            // special clause for DSLDA-NSLT2
            if(phase==0 && option==8 && (j<k1))
                tmpptr[j*ndistWords+i] = 0; // DSLDA-OSST; zero out in test phase
            
        }
        free(tmp1);
        // to be done -- renormalization -- not necessary though
    }
    
    // update sufficient statistics
    for (i = 0; i < nK; i++)
        for (j = 0; j < ndistWords; j++) {
        int index = i + ((int)windexn[j]-1)*nK;    // (j,i)th element
        double value = wcountn[j]*tmpptr[j+i*ndistWords];
        if(index>=nK*V)
        {
            mexPrintf("\n %f %d %d %d ERROR !! \n", windexn[j], i, j, index);
            break;
        }
        sstopicwordptr[index] += value;
        //mexPrintf("%d\t", sstopicwordptr[index], index);
        sstopicptr[i] += value;
        }
    
    mxSetCell (retptr, n, tmp);
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *model, *data;
    mxArray *phi, *phin, *windex, *wcount, *nwordspdoc, *mu, *classlabels, *eta;
    double  *tmp3, *tmp31, *tmp32, *tmp5, *windexn, *wcountn, *annotations, *psigamma;
    double *sstopicwordptr, *sstopicptr;
    int N, nK, C2, Y, V;
    int phase; // 1 for training, 0 for testing
    int option = (int)mxGetScalar(prhs[3]);
    // 1 for unsupervised LDA, 2 for Labeled LDA, 3 for Med-LDA, 4 for DSLDA, 5 for DSLDA with no shared latent topic
    
    model = prhs[0];
    data  = prhs[1];
    //mexPrintf("\nOK....inside update_phi1");
    
    phi = mxGetField(prhs[0],0,"phi");
    nK  = (int)mxGetScalar(mxGetField(prhs[0],0,"K"));
    N   = (int)mxGetScalar(mxGetField(prhs[0],0,"N"));
    V   = (int)mxGetScalar(mxGetField(prhs[0],0,"V"));
    psigamma = mxGetPr(prhs[2]);
    
    option  = (int)mxGetScalar(prhs[3]);
    windex  = mxGetField(prhs[1],0,"windex");
    wcount  = mxGetField(prhs[1],0,"wcount");
    
    if(option>=2)
    {
        phase = (int)mxGetScalar(prhs[4]);
        annotations = mxGetPr(prhs[5]);                        // annotations of visible topics
    }
    
    //if(phase==0)
    //    mexPrintf("\nOK....inside update_phi2");
    plhs[0] = mxCreateCellMatrix(1,N);
    plhs[1] = mxCreateDoubleMatrix(nK,V,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,nK,mxREAL);
    
    sstopicwordptr = mxGetPr(plhs[1]);
    sstopicptr = mxGetPr(plhs[2]);
    
    for (int n = 0; n < N; n++)
    {       
        phin    = mxDuplicateArray (mxGetCell (phi, n));               //pointer to phi{n}
        //if(phase==0)
        //    mexPrintf("\nOK....inside update_phi3");
        int a = mxGetN(mxDuplicateArray (mxGetCell (windex, n)));
        int b = mxGetN(mxDuplicateArray (mxGetCell (wcount, n)));
        //if(phase==0) // debugging for test phase
        //    mexPrintf("doc number %d %d %d\n",n,a,b);
        windexn = mxGetPr(mxDuplicateArray (mxGetCell (windex, n)));   //pointer to windex{n}
        wcountn = mxGetPr(mxDuplicateArray (mxGetCell (wcount, n)));   //pointer to wcount{n}
        update_phi_n(plhs[0], sstopicwordptr, sstopicptr, phin, windexn, wcountn, model, data, psigamma, n, phase, annotations, option);
    }
    
    return;
}
