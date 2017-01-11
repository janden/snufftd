#include "nufftd_spread.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void nufftd_spread(double *tau_re, double *tau_im, int N, int n, int d, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m)
{
    int j, k, l, carry, ind_k, updated;
    int *i, *ind;
    double *P1, *Pj;
    double *Pji;

    int *mu;
    double delta, mult;
    int *mnPowers;

    i = calloc(d, sizeof(int));
    ind = calloc(d+1, sizeof(int));

    P1 = (double *) calloc(q+1, sizeof(double));
    Pj = (double *) calloc(d*(q+1), sizeof(double));

    Pji = (double *) calloc(d+1, sizeof(double));

    mu = calloc(d, sizeof(int));

    mnPowers = calloc(d, sizeof(int));

    for(j = 0; j < q+1; j++)
    {
        P1[j] = exp(-(j-q/2)*(j-q/2)/(4*b));
    }

    mnPowers[0] = 1;
    for(k = 1; k < d; k++)
    {
        mnPowers[k] = mnPowers[k-1]*m*N;
    }

    Pji[d] = 1/pow(2*sqrt(b*M_PI), d);

    updated = d;

    for(j = 0; j < n; j++)
    {
        memset(i, 0, d*sizeof(int));

        for(k = 0; k < d; k++)
        {
            mu[k] = (int) round(m*omega[j+k*n]);
        }

        for(k = 0; k < d; k++)
        {
            delta = m*omega[j+k*n]-mu[k];

            mult = exp(2*delta/(4*b));

            Pj[k+d*(q/2)] = exp(-delta*delta/(4*b));
            Pj[k+d*(q/2+1)] = Pj[k+d*(q/2)]*mult;
            Pj[k+d*(q/2-1)] = Pj[k+d*(q/2)]/mult;

            for(l = 2; l <= q/2; l++)
            {
                Pj[k+d*(q/2+l)] = Pj[k+d*(q/2+l-1)]*mult;
                Pj[k+d*(q/2-l)] = Pj[k+d*(q/2-l+1)]/mult;
            }

            for(l = 0; l < q+1; l++)
            {
                Pj[k+d*l] *= P1[l];
            }
        }

        while(1)
        {
            for(k = updated-1; k > 0; k--)
            {
                ind_k = mu[k]+i[k]-q/2+m*N/2;
                while(ind_k < 0)
                {
                    ind_k += m*N;
                }
                while(ind_k >= m*N)
                {
                    ind_k -= m*N;
                }
                ind[k] = ind_k*mnPowers[k]+ind[k+1];

                Pji[k] = Pji[k+1]*Pj[k+d*i[k]];
            }

            ind_k = mu[0]-q/2+m*N/2;
            while(ind_k < 0)
            {
                ind_k += m*N;
            }
            while(ind_k >= m*N)
            {
                ind_k -= m*N;
            }

            for(i[0] = 0; i[0] < q+1; i[0]++)
            {
                ind[0] = ind_k+ind[1];
                Pji[0] = Pji[1]*Pj[d*i[0]];

                tau_re[ind[0]] += Pji[0]*alpha_re[j];
                if(alpha_im != NULL)
                {
                    tau_im[ind[0]] += Pji[0]*alpha_im[j];
                }

                ind_k++;
                if(ind_k == m*N)
                {
                    ind_k = 0;
                }
            }

            carry = 1;
            k = 1;
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
    free(Pj);

    free(Pji);

    free(mu);

    free(mnPowers);
}
