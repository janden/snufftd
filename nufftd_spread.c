#include "nufftd_spread.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void nufftd_spread(double *tau_re, double *tau_im, int N, int n, int d, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m)
{
    int j, k, l, carry, ind_k, updated;
    int *i, *ind;
    double *P1, *P2, *P3;
    double *Pji;

    int *mu;
    double Pji0, delta;
    int *mnPowers;

    i = calloc(d, sizeof(int));
    ind = calloc(d+1, sizeof(int));

    P1 = (double *) calloc(q+1, sizeof(double));
    P2 = (double *) calloc(d, sizeof(double));
    P3 = (double *) calloc(d, sizeof(double));

    Pji = (double *) calloc(d+1, sizeof(double));

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

    Pji[d] = Pji0;

    updated = d;

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
            for(k = updated-1; k >= 0; k--)
            {
                ind_k = (mu[j+k*n]+i[k]-q/2+m*N/2) % (m*N);
                ind_k = (ind_k < 0) ? (ind_k+m*N) : ind_k;
                ind[k] = ind_k*mnPowers[k]+ind[k+1];

                Pji[k] = Pji[k+1]*P1[i[k]]*P2[k];
                if(i[k] < q/2)
                {
                    for(l = 0; l < q/2-i[k]; l++)
                    {
                        Pji[k] /= P3[k];
                    }
                }
                else
                {
                    for(l = 0; l < i[k]-q/2; l++)
                    {
                        Pji[k] *= P3[k];
                    }
                }
            }

            tau_re[ind[0]] += Pji[0]*alpha_re[j];
            if(alpha_im != NULL)
            {
                tau_im[ind[0]] += Pji[0]*alpha_im[j];
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
            updated = k;

            if(carry)
            {
                break;
            }
        }
    }

    free(i);
    free(ind);

    free(P1);
    free(P2);
    free(P3);

    free(Pji);

    free(mu);

    free(mnPowers);
}
