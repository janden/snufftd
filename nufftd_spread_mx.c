#include <mex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

void checkInputsOutputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs != 6)
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

    if(!mxIsNumeric(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1)
    {
        mexErrMsgTxt("b must be a numeric scalar.");
    }

    if(mxGetScalar(prhs[3]) < 0)
    {
        mexErrMsgTxt("b must be positive.");
    }

    if(!mxIsNumeric(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1)
    {
        mexErrMsgTxt("q must be a numeric scalar.");
    }

    if(fabs(mxGetScalar(prhs[4])-round(mxGetScalar(prhs[4]))) > 1e-15 || mxGetScalar(prhs[4]) < 1)
    {
        mexErrMsgTxt("q must be a positive integer.");
    }

    if(!mxIsNumeric(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1)
    {
        mexErrMsgTxt("m must be a numeric scalar.");
    }

    if(fabs(mxGetScalar(prhs[5])-round(mxGetScalar(prhs[5]))) > 1e-15 || mxGetScalar(prhs[5]) < 1)
    {
        mexErrMsgTxt("m must be a positive integer.");
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

void parseInputs(const mxArray *prhs[], mwSize *N, mwSize *n, mwSize *d, double **omega, double **alpha_re, double **alpha_im, double *b, int *q, int *m)
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

    *b = mxGetScalar(prhs[3]);
    *q = (int) round(mxGetScalar(prhs[4]));
    *m = (int) round(mxGetScalar(prhs[5]));
}

void prepareOutputs(mxArray *plhs[], mwSize N, mwSize d, int m, double **tau_re, double **tau_im)
{
    mwSize *dims;
    int i, ndims;

    ndims = (d > 1 ? d : 2);

    dims = calloc(ndims, sizeof(mwSize));

    for(i = 0; i < ndims; i++)
    {
        dims[i] = (i < d ? m*N : 1);
    }

    plhs[0] = mxCreateNumericArray(d > 1 ? d : 2, dims, mxDOUBLE_CLASS, mxCOMPLEX);

    free(dims);

    *tau_re = mxGetPr(plhs[0]);
    *tau_im = mxGetPi(plhs[0]);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize N, n, d;
    double *omega, *alpha_re, *alpha_im, *tau_re, *tau_im;

    double b;
    int q, m;

    checkInputsOutputs(nlhs, plhs, nrhs, prhs);

    parseInputs(prhs, &N, &n, &d, &omega, &alpha_re, &alpha_im, &b, &q, &m);

    prepareOutputs(plhs, N, d, m, &tau_re, &tau_im);
}
