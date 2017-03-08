#include <mex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "nufftd_spread.h"

void check_inputs_outputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i;

    if(nrhs != 7)
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

    if(!mxIsDouble(prhs[1]) || mxIsSparse(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) > 2 || mxIsComplex(prhs[1]) ||
       mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("grid_shift must be a full double matrix of size d-by-1.");
    }

    if(!mxIsDouble(prhs[2]) || mxIsSparse(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) > 2 || mxIsComplex(prhs[2]))
    {
        mexErrMsgTxt("omega must be a full real double matrix of size n-by-d.");
    }

    if(!mxIsDouble(prhs[3]) || mxIsSparse(prhs[3]) || mxGetNumberOfDimensions(prhs[3]) > 2 || mxGetN(prhs[3]) != 1)
    {
        mexErrMsgTxt("alpha must be a full double matrix of size n-by-1.");
    }

    if(!mxIsNumeric(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1)
    {
        mexErrMsgTxt("b must be a numeric scalar.");
    }

    if(mxGetScalar(prhs[4]) < 0)
    {
        mexErrMsgTxt("b must be positive.");
    }

    if(!mxIsNumeric(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1)
    {
        mexErrMsgTxt("q must be a numeric scalar.");
    }

    if(fabs(mxGetScalar(prhs[5])-round(mxGetScalar(prhs[5]))) > 1e-15 || mxGetScalar(prhs[5]) < 1)
    {
        mexErrMsgTxt("q must be a positive integer.");
    }

    if(!mxIsNumeric(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1)
    {
        mexErrMsgTxt("m must be a numeric scalar.");
    }

    if(fabs(mxGetScalar(prhs[6])-round(mxGetScalar(prhs[6]))) > 1e-15 || mxGetScalar(prhs[6]) < 1)
    {
        mexErrMsgTxt("m must be a positive integer.");
    }

    if(mxGetM(prhs[1]) != mxGetN(prhs[2]))
    {
        mexErrMsgTxt("size(grid_shift, 1) must match size(omega, 2).");
    }

    if(mxGetM(prhs[2]) != mxGetM(prhs[3]))
    {
        mexErrMsgTxt("size(omega, 1) must match size(alpha, 1).");
    }

    for(i = 0; i < mxGetM(prhs[1]); i++)
    {
        if(fabs(mxGetPr(prhs[1])[i]-round(mxGetPr(prhs[1])[i])) > 1e-15 || mxGetPr(prhs[1])[i] < 0 ||
           mxGetPr(prhs[1])[i] >= mxGetScalar(prhs[6]))
        {
            mexErrMsgTxt("grid_shift must be composed of integers between 0 and m-1.");
        }
    }

    if(nlhs > 1)
    {
        mexErrMsgTxt("Incorrect number of outputs.");
    }
}

void parse_inputs(const mxArray *prhs[], mwSize *N, mwSize *n, mwSize *d, int **grid_shift, double **omega, double **alpha_re, double **alpha_im, double *b, int *q, int *m)
{
    int i;

    *N = (mwSize) round(mxGetScalar(prhs[0]));

    *d = mxGetM(prhs[1]);

    *grid_shift = calloc(*d, sizeof(int));

    for(i = 0; i < *d; i++)
    {
        (*grid_shift)[i] = (int) round(mxGetPr(prhs[1])[i]);
    }

    *n = mxGetM(prhs[2]);

    *omega = mxGetPr(prhs[2]);

    *alpha_re = mxGetPr(prhs[3]);

    if(mxIsComplex(prhs[3]))
    {
        *alpha_im = mxGetPi(prhs[3]);
    }
    else
    {
        *alpha_im = NULL;
    }

    *b = mxGetScalar(prhs[4]);
    *q = (int) round(mxGetScalar(prhs[5]));
    *m = (int) round(mxGetScalar(prhs[6]));
}

void prepare_outputs(mxArray *plhs[], mwSize N, mwSize d, int m, double **tau_re, double **tau_im)
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

    *tau_re = mxGetPr(plhs[0]);
    *tau_im = mxGetPi(plhs[0]);
}

void free_inputs(int *grid_shift)
{
    free(grid_shift);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize N, n, d;
    double *omega, *alpha_re, *alpha_im, *tau_re, *tau_im;
    int *grid_shift;

    double b;
    int q, m;

    int i;

    check_inputs_outputs(nlhs, plhs, nrhs, prhs);

    parse_inputs(prhs, &N, &n, &d, &grid_shift, &omega, &alpha_re, &alpha_im, &b, &q, &m);

    prepare_outputs(plhs, N, d, m, &tau_re, &tau_im);

    sub_snufftd_spread(tau_re, tau_im, N, n, d, grid_shift, omega, alpha_re, alpha_im, b, q, m, NULL);

    free_inputs(grid_shift);
}
