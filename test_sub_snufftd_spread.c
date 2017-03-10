#include <stdlib.h>

#include "nufftd_spread.h"

#define         N       64
#define         D       2
#define         M       16*N*N

int main()
{
    double *tau_re, *tau_im, *omega, *alpha_re, *alpha_im;

    int *grid_shift;

    int k, i;

    if(D == 1)
    {
        tau_re = calloc(N, sizeof(double));
        tau_im = calloc(N, sizeof(double));
    }
    else if(D == 2)
    {
        tau_re = calloc(N*N, sizeof(double));
        tau_im = calloc(N*N, sizeof(double));
    }
    else if(D == 3)
    {
        tau_re = calloc(N*N*N, sizeof(double));
        tau_im = calloc(N*N*N, sizeof(double));
    }

    omega = calloc(M*D, sizeof(double));

    alpha_re = calloc(M, sizeof(double));
    alpha_im = calloc(M, sizeof(double));

    grid_shift = calloc(D, sizeof(int));

    srand(0);

    for(i = 0; i < D; i++)
    {
        grid_shift[i] = 1;
    }

    for(i = 0; i < 400; i++)
    {
        for(k = 0; k < D*M; k++)
        {
            omega[k] = N/2*(((double) rand())/RAND_MAX-0.5);
        }

        for(k = 0; k < M; k++)
        {
            alpha_re[k] = ((double) rand())/RAND_MAX;
            alpha_im[k] = ((double) rand())/RAND_MAX;
        }

        sub_snufftd_spread(tau_re, tau_im, N, M, D, grid_shift, omega, alpha_re, alpha_im, 1.5629, 28, 2, NULL);
    }

    free(tau_re);
    free(tau_im);
    free(omega);
    free(alpha_re);
    free(alpha_im);

    return 0;
}
