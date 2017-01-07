#include <stdlib.h>
#include <math.h>

#include "nudftd.h"

void nudftd(double *f_re, double *f_im, int N, int n, int d, double *omega, double *alpha_re, double *alpha_im)
{
    int j, k, carry;

    int *i;

    size_t ind;

    i = calloc(d, sizeof(int));

    ind = 0;

    while(1)
    {
        for(j = 0; j < n; j++)
        {
            double gamma, co, si;

            gamma = 0;
            for(k = 0; k < d; k++)
            {
                gamma += (i[k]-N/2)*omega[j+k*n];
            }
            gamma = 2*M_PI*gamma/N;

            co = cos(gamma);
            si = sin(gamma);

            if(alpha_im == NULL)
            {
                f_re[ind] += co*alpha_re[j];
                f_im[ind] += si*alpha_re[j];
            }
            else
            {
                f_re[ind] += co*alpha_re[j]-si*alpha_im[j];
                f_im[ind] += si*alpha_re[j]+co*alpha_im[j];
            }
        }

        carry = 1;
        k = 0;
        while(carry > 0 && k < d)
        {
            i[k]++;
            carry = 0;

            if(i[k] == N)
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

        ind++;
    }

    free(i);
}
