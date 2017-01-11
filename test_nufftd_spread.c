#include <stdlib.h>

#include "nufftd_spread.h"

#define         N       32
#define         D       3
#define         M       1024

int main()
{
    double tau_re[(2*N)*(2*N)*(2*N)], tau_im[(2*N)*(2*N)*(2*N)];

    double omega[D*M];

    double alpha_re[M], alpha_im[M];

    int k;

    srand(0);

    for(k = 0; k < D*M; k++)
    {
        omega[k] = N/2*(rand()/RAND_MAX-0.5);
    }

    for(k = 0; k < M; k++)
    {
        alpha_re[k] = rand()/RAND_MAX;
        alpha_im[k] = rand()/RAND_MAX;
    }

    nufftd_spread(tau_re, tau_im, N, M, D, omega, alpha_re, alpha_im, 1.5629, 28, 2);

    return 0;
}
