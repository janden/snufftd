#include "nufftd_spread.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void nufftd_spread(double *tau_re, double *tau_im, int N, int n, int d, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m)
{
    int j, k, l, carry, ind, ind_k;
    int *i;
    double *P1, *P2, *P3;
    double Pji;

    int *mu;
    double Pji0, delta;
    int *mnPowers;

    i = calloc(d, sizeof(int));

    P1 = (double *) calloc(q+1, sizeof(double));
    P2 = (double *) calloc(d, sizeof(double));
    P3 = (double *) calloc(d, sizeof(double));

    mu = calloc(n*d, sizeof(int));

    mnPowers = calloc(d, sizeof(int));

    for(j = 0; j < q+1; j++)
    {
        P1[j] = exp(-(j-q/2)*(j-q/2)/(4*b));
    }

    for(j = 0; j < n*d; j++)
    {
        mu[j] = (int) round(m*omega[j]);
    }

    Pji0 = 1/pow(2*sqrt(b*M_PI), d);

    mnPowers[0] = 1;
    for(k = 1; k < d; k++)
    {
        mnPowers[k] = mnPowers[k-1]*m*N;
    }

    for(j = 0; j < n; j++)
    {
        memset(i, 0, d*sizeof(int));

        for(k = 0; k < d; k++)
        {
            delta = m*omega[j+k*n]-mu[j+k*n];

            P2[k] = exp(-delta*delta/(4*b));
            P3[k] = exp(2*delta/(4*b));
        }

        while(1)
        {
            ind = 0;
            Pji = Pji0;

            for(k = 0; k < d; k++)
            {
                ind_k = (mu[j+k*n]+i[k]-q/2+m*N/2) % (m*N);
                ind_k = (ind_k < 0) ? (ind_k+m*N) : ind_k;
                ind += ind_k*mnPowers[k];

                Pji *= P1[i[k]]*P2[k]*pow(P3[k], i[k]-q/2);
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

    free(mu);

    free(mnPowers);
}
