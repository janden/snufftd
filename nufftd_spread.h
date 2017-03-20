#ifndef NUFFTD_SPREAD_H
#define NUFFTD_SPREAD_H

#include <stdlib.h>

void nufftd_spread(double *tau_re, double *tau_im, size_t N, size_t n, size_t d, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m);

void sub_snufftd_spread(double *tau_re, double *tau_im, size_t N, size_t n, size_t d, int *grid_shift, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m, double *precomp);

#endif
