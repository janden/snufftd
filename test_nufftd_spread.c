#include <stdlib.h>

#include "nufftd_spread.h"

#define         N       64
#define         D       2
#define         M       16*N*N

int main()
{
    double *tau_re, *tau_im, *omega, *alpha_re, *alpha_im;

    int k, i;

    if(D == 1)
    {
        tau_re = calloc(2*N, sizeof(double));
        tau_im = calloc(2*N, sizeof(double));
    }
    else if(D == 2)
    {
        tau_re = calloc(2*N*2*N, sizeof(double));
        tau_im = calloc(2*N*2*N, sizeof(double));
    }
    else if(D == 3)
    {
        tau_re = calloc(2*N*2*N*2*N, sizeof(double));
        tau_im = calloc(2*N*2*N*2*N, sizeof(double));
    }

    omega = calloc(M*D, sizeof(double));

    alpha_re = calloc(M, sizeof(double));
    alpha_im = calloc(M, sizeof(double));

    srand(0);

    for(i = 0; i < 100; i++)
    {
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
    }

    free(tau_re);
    free(tau_im);
    free(omega);
    free(alpha_re);
    free(alpha_im);

    return 0;
}
