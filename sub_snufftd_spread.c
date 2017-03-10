#include "nufftd_spread.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void sub_snufftd_spread(double *tau_re, double *tau_im, int N, int n, int d, int *grid_shift, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m, double *precomp)
{
    int k, l, carry, ind_k, updated;
    size_t j;
    int *i;
    size_t *ind;
    double *P1, *Pj;
    double *Pji;

    int *mu;
    double delta, mult, mult1, val, alpha_re_j, alpha_im_j;
    size_t *nPowers;

    int max_width, center_ind;
    int *mu_shift;

    int *num_ind;

    int extm, extp;

    max_width = 2*(q/(2*m))+2;
    center_ind = q/(2*m);

    i = calloc(d, sizeof(int));
    ind = calloc(d+1, sizeof(size_t));

    P1 = (double *) calloc(q+1, sizeof(double));
    Pj = (double *) calloc(d*max_width, sizeof(double));

    Pji = (double *) calloc(d+1, sizeof(double));

    mu = calloc(d, sizeof(int));

    nPowers = calloc(d, sizeof(size_t));

    mu_shift = calloc(d, sizeof(int));

    num_ind = calloc(d, sizeof(int));

    for(l = 0; l < q+1; l++)
    {
        P1[l] = exp(-(l-q/2)*(l-q/2)/(4*b));
    }

    nPowers[0] = (size_t) 1;
    for(k = 1; k < d; k++)
    {
        nPowers[k] = nPowers[k-1]*((size_t) N);
    }

    Pji[d] = 1/pow(2*sqrt(b*M_PI), d);

    updated = d;

    for(j = 0; j < n; j++)
    {
        for(k = 0; k < d; k++)
        {
            i[k] = 0;
        }

        for(k = 0; k < d; k++)
        {
            mu[k] = (int) round(m*omega[k+j*d]);
        }

        for(k = 0; k < d; k++)
        {
            mu_shift[k] = (mu[k]-grid_shift[k])%m;
            if(mu_shift[k] < 0)
            {
                mu_shift[k] = mu_shift[k]+m;
            }
        }

        for(k = 0; k < d; k++)
        {
            extm = (q/2-mu_shift[k])/m;
            extp = (q/2+mu_shift[k])/m;

            delta = m*omega[k+j*d]-mu[k];

            if(precomp == NULL)
            {
                val = exp(-delta*delta/(4*b));

                mult1 = exp(-2*m*delta/(4*b));
            }
            else
            {
                val = precomp[0+2*k+2*d*j];

                mult1 = precomp[1+2*k+2*d*j];
            }

            for(l = 0; l < mu_shift[k]; l++)
            {
                val = val*mult1;
            }

            mult = mult1;

            for(l = 1; l < m; l++)
            {
                mult = mult*mult1;
            }

            Pj[extm+k*max_width] = val;

            for(l = 1; l <= extm; l++)
            {
                val = val*mult;
                Pj[extm-l+k*max_width] = P1[q/2-mu_shift[k]-m*l]*val;
            }

            val = Pj[extm+k*max_width];
            mult = 1/mult;

            for(l = 1; l <= extp; l++)
            {
                val = val*mult;
                Pj[extm+l+k*max_width] = P1[q/2-mu_shift[k]+m*l]*val;
            }

            Pj[extm+k*max_width] = P1[q/2-mu_shift[k]]*Pj[extm+k*max_width];

            mu[k] = (mu[k]-grid_shift[k]-mu_shift[k])/m-extm+N/2;

            while(mu[k] < 0)
            {
                mu[k] += N;
            }

            num_ind[k] = extm+extp+1;
        }

        while(1)
        {
            for(k = updated-1; k > 0; k--)
            {
                ind_k = mu[k]+i[k];
                while(ind_k >= N)
                {
                    ind_k -= N;
                }
                ind[k] = ((size_t) ind_k)*nPowers[k]+ind[k+1];

                Pji[k] = Pji[k+1]*Pj[i[k]+max_width*k];
            }

            alpha_re_j = Pji[1]*alpha_re[j];
            if(alpha_im != NULL)
            {
                alpha_im_j = Pji[1]*alpha_im[j];
            }

            ind[0] = ((size_t) mu[k])+ind[1];

            for(i[0] = 0; i[0] < num_ind[0]; i[0]++)
            {
                tau_re[ind[0]] += Pj[i[0]]*alpha_re_j;
                if(alpha_im != NULL)
                {
                    tau_im[ind[0]] += Pj[i[0]]*alpha_im_j;
                }

                ind[0]++;
                if(ind[0] == ind[1]+N)
                {
                    ind[0] = ind[1];
                }
            }

            carry = 1;
            k = 1;
            while(carry > 0 && k < d)
            {
                i[k]++;
                carry = 0;

                if(i[k] == num_ind[k])
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

    free(nPowers);

    free(mu_shift);

    free(num_ind);
}

