#include "nufftd_spread.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void nufftd_spread(double *tau_re, double *tau_im, size_t N, size_t n, size_t d, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m)
{
    int k, l, carry, ind_k, updated, stop;
    size_t j;
    int *i;
    size_t *ind;
    double *P1, *Pj;
    double *Pji;

    int *mu;
    double delta, mult, val, alpha_re_j, alpha_im_j;
    size_t *mnPowers;

    i = calloc(d, sizeof(int));
    ind = calloc(d+1, sizeof(size_t));

    P1 = (double *) calloc(q+1, sizeof(double));
    Pj = (double *) calloc(d*(q+1), sizeof(double));

    Pji = (double *) calloc(d+1, sizeof(double));

    mu = calloc(d, sizeof(int));

    mnPowers = calloc(d, sizeof(size_t));

    for(l = 0; l < q+1; l++)
    {
        P1[l] = exp(-(l-q/2)*(l-q/2)/(4*b));
    }

    mnPowers[0] = (size_t) 1;
    for(k = 1; k < d; k++)
    {
        mnPowers[k] = mnPowers[k-1]*((size_t) (m*N));
    }

    Pji[d] = 1/pow(2*sqrt(b*M_PI), d);

    updated = d;

    for(j = 0; j < n; j++)
    {
        for(k = 0; k < d; k++)
        {
            mu[k] = (int) round(m*omega[k+j*d]);
        }

        for(k = 0; k < d; k++)
        {
            delta = m*omega[k+j*d]-mu[k];

            val = exp(-delta*delta/(4*b));

            mult = exp(2*delta/(4*b));

            Pj[(q/2)+k*(q+1)] = val;

            for(l = 1; l <= q/2; l++)
            {
                val = val*mult;
                Pj[(q/2+l)+k*(q+1)] = P1[q/2+l]*val;
            }

            val = Pj[(q/2)+k*(q+1)];
            mult = 1/mult;

            for(l = 1; l <= q/2; l++)
            {
                val = val*mult;
                Pj[(q/2-l)+k*(q+1)] = P1[q/2-l]*val;
            }

            mu[k] = mu[k]-q/2+m*N/2;

            while(mu[k] < 0)
            {
                mu[k] += m*N;
            }
        }

        memset(i, 0, d*sizeof(int));

        while(1)
        {
            for(k = updated-1; k > 0; k--)
            {
                ind_k = mu[k]+i[k];
                while(ind_k >= m*N)
                {
                    ind_k -= m*N;
                }
                ind[k] = ((size_t) ind_k)*mnPowers[k]+ind[k+1];

                Pji[k] = Pji[k+1]*Pj[i[k]+(q+1)*k];
            }

            alpha_re_j = Pji[1]*alpha_re[j];
            if(alpha_im != NULL)
            {
                alpha_im_j = Pji[1]*alpha_im[j];
            }

            ind[0] = ((size_t) mu[0])+ind[1];

            i[0] = 0;

            while(1)
            {
                stop = q+1;

                if(ind[0]+stop-i[0] > ind[1]+m*N)
                {
                    stop = ind[1]+m*N-ind[0]+i[0];
                }

                while(i[0] < stop)
                {
                    tau_re[ind[0]] += Pj[i[0]]*alpha_re_j;
                    if(alpha_im != NULL)
                    {
                        tau_im[ind[0]] += Pj[i[0]]*alpha_im_j;
                    }

                    ind[0]++;
                    i[0]++;
                }

                if(i[0] == q+1)
                {
                    break;
                }

                ind[0] = ind[1];
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
