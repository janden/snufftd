#ifndef NUFFTD_SPREAD_H
#define NUFFTD_SPREAD_H

void nufftd_spread(double *tau_re, double *tau_im, int N, int n, int d, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m);

void sub_snufftd_spread(double *tau_re, double *tau_im, int N, int n, int d, int *grid_shift, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m, double *precomp);

#endif
