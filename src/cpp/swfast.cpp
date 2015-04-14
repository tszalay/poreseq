/*==========================================================
 *
 * swfast.cpp
 *
 * Algorithm for computing a striped Smith-Waterman alignment
 * when the rough alignment is known through other means
 * (or just for computing using full matrix otherwise)
 *
 * Arguments: seq0, seq1, paired inds, stripe width
 *
 *========================================================*/

#include "matrix.h"
#include "mex.h"

#include "swlib.h"

using namespace std;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2)
        mexErrMsgIdAndTxt("CrampBio:swfast","Not enough arguments.");

    // read in sequences
    if (!mxIsChar(prhs[0]) || !mxIsChar(prhs[1]))
        mexErrMsgIdAndTxt("CrampBio:swfast","Arrays not of type char.");
    
    int n0 = mxGetNumberOfElements(prhs[0]);
    int n1 = mxGetNumberOfElements(prhs[1]);
    
    uint16_t* seq0 = (uint16_t*)mxGetPr(prhs[0]);
    uint16_t* seq1 = (uint16_t*)mxGetPr(prhs[1]);
    
    SWAlignment align;

    // did we receive striped indices? if not, do full matrix
    if (nrhs > 2)
    {
        int stripe_width = (int)mxGetScalar(prhs[3]);
        
        double* pr = (double*)mxGetPr(prhs[2]);
        // calculate linear stripe coefficients from given points
        // the order is pr = [i0 j0 i1 j1]
        double align_m = (pr[2]-pr[0])/(double)(pr[3]-pr[1]);
        double align_b = pr[0] - align_m*pr[1];
        
        align = swfast(vector<uint8_t>(seq0,seq0+n0),vector<uint8_t>(seq1,seq1+n1),
                align_m, align_b, stripe_width);
    }
    else
    {
        align = swfull(vector<uint8_t>(seq0,seq0+n0),vector<uint8_t>(seq1,seq1+n1));
    }

    // save to mxarray
    plhs[0] = mxCreateDoubleScalar(align.score);
    if (nlhs > 1)
    {
        int n = align.inds1.size();
        plhs[1] = mxCreateDoubleMatrix(n,2,mxREAL);
        double *pr = mxGetPr(plhs[1]);

        for (int i=0; i<n; i++)
        {
            pr[i] = align.inds1[i];
            pr[i+n] = align.inds2[i];
        }
    }
}



