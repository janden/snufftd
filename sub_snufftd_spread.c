#include "nufftd_spread.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void sub_snufftd_spread(double *tau_re, double *tau_im, int N, int n, int d, int *grid_shift, double *omega, double *alpha_re, double *alpha_im, double b, int q, int m)
{
    int j, k, l, carry, ind_k, updated;
    int *i, *ind;
    double *P1, *Pj;
    double *Pji;

    int *mu;
    double delta, mult, val, alpha_re_j, alpha_im_j;
    int *nPowers;

    int max_width, center_ind;
    int *mu_shift;

    int *start_ind, *num_ind;

    max_width = 2*(q/(2*m))+2;
    center_ind = q/(2*m);

    i = calloc(d, sizeof(int));
    ind = calloc(d+1, sizeof(int));

    P1 = (double *) calloc(q+1, sizeof(double));
    Pj = (double *) calloc(d*max_width, sizeof(double));

    Pji = (double *) calloc(d+1, sizeof(double));

    mu = calloc(d, sizeof(int));

    nPowers = calloc(d, sizeof(int));

    mu_shift = calloc(d, sizeof(int));

    start_ind = calloc(d, sizeof(int));
    num_ind = calloc(d, sizeof(int));

    for(j = 0; j < q+1; j++)
    {
        P1[j] = exp(-(j-q/2)*(j-q/2)/(4*b));
    }

    nPowers[0] = 1;
    for(k = 1; k < d; k++)
    {
        nPowers[k] = nPowers[k-1]*N;
    }

    Pji[d] = 1/pow(2*sqrt(b*M_PI), d);

    updated = d;

    for(j = 0; j < n; j++)
    {
        memset(i, 0, d*sizeof(int));

        for(k = 0; k < d; k++)
        {
            mu[k] = (int) round(m*omega[j+k*n]);
            mu_shift[k] = (mu[k]-grid_shift[k])%m;
            if(mu_shift[k] < 0)
            {
                mu_shift[k] = mu_shift[k]+m;
            }
        }

        for(k = 0; k < d; k++)
        {
            delta = m*omega[j+k*n]-mu[k];

            val = exp(-(delta*delta+2*mu_shift[k]*delta)/(4*b));

            Pj[center_ind+k*max_width] = val;

            mult = exp(-2*m*delta/(4*b));

            for(l = 1; l <= (q/2-mu_shift[k])/m; l++)
            {
                val = val*mult;
                Pj[center_ind-l+k*max_width] = P1[q/2-mu_shift[k]-m*l]*val;
            }

            val = Pj[center_ind+k*max_width];
            mult = 1/mult;

            for(l = 1; l <= (q/2+mu_shift[k])/m; l++)
            {
                val = val*mult;
                Pj[center_ind+l+k*max_width] = P1[q/2-mu_shift[k]+m*l]*val;
            }

            Pj[center_ind+k*max_width] = P1[q/2-mu_shift[k]]*Pj[center_ind+k*max_width];

            start_ind[k] = -(q/2-mu_shift[k])/m;
            num_ind[k] = (q/2-mu_shift[k])/m+(q/2+mu_shift[k])/m+1;
        }

        while(1)
        {
            for(k = updated-1; k > 0; k--)
            {
                ind_k = (mu[k]-grid_shift[k]-mu_shift[k])/m+i[k]+start_ind[k]+N/2;
                while(ind_k < 0)
                {
                    ind_k += N;
                }
                while(ind_k >= N)
                {
                    ind_k -= N;
                }
                ind[k] = ind_k*nPowers[k]+ind[k+1];

                Pji[k] = Pji[k+1]*Pj[i[k]+start_ind[k]+center_ind+max_width*k];
            }

            alpha_re_j = Pji[1]*alpha_re[j];
            if(alpha_im != NULL)
            {
                alpha_im_j = Pji[1]*alpha_im[j];
            }

            ind_k = (mu[0]-grid_shift[0]-mu_shift[0])/m+start_ind[0]+N/2;

            while(ind_k < 0)
            {
                ind_k += N;
            }
            while(ind_k >= N)
            {
                ind_k -= N;
            }

            for(i[0] = 0; i[0] < num_ind[0]; i[0]++)
            {
                ind[0] = ind_k+ind[1];

                tau_re[ind[0]] += Pj[i[0]+start_ind[0]+center_ind]*alpha_re_j;
                if(alpha_im != NULL)
                {
                    tau_im[ind[0]] += Pj[i[0]+start_ind[0]+center_ind]*alpha_im_j;
                }

                ind_k++;
                if(ind_k == N)
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

    free(start_ind);
    free(num_ind);
}

