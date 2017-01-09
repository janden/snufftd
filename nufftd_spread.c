#include "nufftd_spread.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void precalc_gaussian_kernel(int n, int d, double *omega, double b, int q, int m, double *P1, double *P2, double *P3)
{
    int j;
    double delta;

    for(j = 0; j < q+1; j++)
    {
        P1[j] = exp(-(j-q/2)*(j-q/2)/(4*b));
    }

    for(j = 0; j < n*d; j++)
    {
        delta = m*omega[j]-round(m*omega[j]);

        P2[j] = exp(-delta*delta/(4*b));
        P3[j] = exp(2*delta/(4*b));
    }
}

void nufftd_spread(double *tau_re, double *tau_im, int N, int n, int d, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m)
{
    int j, k, carry, ind, ind_k;
    int *i;
    double *P1, *P2, *P3;
    double Pji;

    i = calloc(d, sizeof(int));

    P1 = (double *) calloc(q+1, sizeof(double));
    P2 = (double *) calloc(n*d, sizeof(double));
    P3 = (double *) calloc(n*d, sizeof(double));

    precalc_gaussian_kernel(n, d, omega, b, q, m, P1, P2, P3);

    for(j = 0; j < n; j++)
    {
        memset(i, 0, d*sizeof(int));

        while(1)
        {
            ind = 0;
            Pji = 1/pow(2*sqrt(b*M_PI), d);

            for(k = 0; k < d; k++)
            {
                ind_k = (((int) round(m*omega[j+k*n]))+i[k]-q/2+m*N/2) % (m*N);
                ind_k = (ind_k < 0) ? (ind_k+m*N) : ind_k;
                ind += ind_k*pow(m*N, k);

                Pji *= P1[i[k]]*P2[j+k*n]*pow(P3[j+k*n], i[k]-q/2);
            }

            tau_re[ind] += Pji*alpha_re[j];
            if(alpha_im != NULL)
            {
                tau_im[ind] += Pji*alpha_im[j];
            }

            carry = 1;
            k = 0;
            while(carry > 0 && k < d)
            {
                i[k]++;
                carry = 0;

                if(i[k] == q+1)
                {
                    i[k] = 0;
                    carry = 1;
                }

                k++;
            }

            if(carry)
            {
                break;
            }
        }
    }

    free(i);

    free(P1);
    free(P2);
    free(P3);
}
