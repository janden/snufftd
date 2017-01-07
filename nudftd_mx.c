#include <mex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "nudftd.h"

void checkInputsOutputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs != 3)
    {
        mexErrMsgTxt("Incorrect number of inputs.");
    }

    if(!mxIsNumeric(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 1)
    {
        mexErrMsgTxt("N must be a numeric scalar.");
    }

    if(fabs(mxGetScalar(prhs[0])-round(mxGetScalar(prhs[0]))) > 1e-15 || mxGetScalar(prhs[0]) < 1)
    {
        mexErrMsgTxt("N must be a positive integer.");
    }

    if(!mxIsDouble(prhs[1]) || mxIsSparse(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) > 2 || mxIsComplex(prhs[1]))
    {
        mexErrMsgTxt("omega must be a full real double matrix of size n-by-d.");
    }

    if(!mxIsDouble(prhs[2]) || mxIsSparse(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) > 2 || mxGetN(prhs[2]) != 1)
    {
        mexErrMsgTxt("alpha must be a full double matrix of size n-by-1.");
    }

    if(mxGetM(prhs[1]) != mxGetM(prhs[2]))
    {
        mexErrMsgTxt("size(omega, 1) must match size(alpha, 1).");
    }

    if(nlhs > 1)
    {
        mexErrMsgTxt("Incorrect number of outputs.");
    }
}

void parseInputs(const mxArray *prhs[], mwSize *N, mwSize *n, mwSize *d, double **omega, double **alpha_re, double **alpha_im)
{
    *N = (mwSize) round(mxGetScalar(prhs[0]));
    *n = mxGetM(prhs[1]);
    *d = mxGetN(prhs[1]);

    *omega = mxGetPr(prhs[1]);

    *alpha_re = mxGetPr(prhs[2]);

    if(mxIsComplex(prhs[2]))
    {
        *alpha_im = mxGetPi(prhs[2]);
    }
    else
    {
        *alpha_im = NULL;
    }
}

void prepareOutputs(mxArray *plhs[], mwSize N, mwSize d, double **f_re, double **f_im)
{
    mwSize *dims;
    int i, ndims;

    ndims = (d > 1 ? d : 2);

    dims = calloc(ndims, sizeof(mwSize));

    for(i = 0; i < ndims; i++)
    {
        dims[i] = (i < d ? N : 1);
    }

    plhs[0] = mxCreateNumericArray(d > 1 ? d : 2, dims, mxDOUBLE_CLASS, mxCOMPLEX);

    free(dims);

    *f_re = mxGetPr(plhs[0]);
    *f_im = mxGetPi(plhs[0]);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize N, n, d;
    double *omega, *alpha_re, *alpha_im, *f_re, *f_im;

    checkInputsOutputs(nlhs, plhs, nrhs, prhs);

    parseInputs(prhs, &N, &n, &d, &omega, &alpha_re, &alpha_im);

    prepareOutputs(plhs, N, d, &f_re, &f_im);

    nudftd(f_re, f_im, N, n, d, omega, alpha_re, alpha_im);
}
