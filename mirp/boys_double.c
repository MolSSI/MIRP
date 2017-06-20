/*! \file
 *
 * \brief Calculation of the boys function in double precision
 */

#include "mirp/boys.h"
#include "mirp/math.h"
#include <math.h>
#include <assert.h>

static void mirp_boys_double_nonzero(double *F, int m, double t)
{
    int i;
    double test;
    double term;
    double sum;
    double oldterm;

    const double t2 = 2 * t;
    const double et = exp(-t);

    int do_short = 0;


    /* The short-range formula converges much better for
     * t < (m + 3/2)
     * So skip the long range if that happens
     * Note that this is a conservative bound, and could probably be increased
     */
    if(t < ((double)m + 0.5))
        do_short = 1;
     

    if(!do_short)
    {
        /* Attempt the long-range approximation */
        F[m] = mirp_double_factorial(2*m-1) / pow(2, m+1);
        F[m] *= sqrt(ARBINT_PI / pow(t, 2*m+1));

        /* Determine the error associated with the long-range approximation */
        term = 1.0;
        sum = 0.0;

        i = 0;
        do {
            i++;
            oldterm = fabs(term);
            term *= (2*m - 2*i + 1) / t2;

            test = fabs(term);
            if(test > oldterm)
               break; 

            test = sum;
            sum += term;
            test -= sum;
        } while(test != 0.0);

        sum *= et / t2;

        /*
         * Determine if this error is satisfactory
         * If not, mark that we have to do the short-range version
         */
        test = F[m] - sum;
        if( (test - F[m]) != 0.0 )
            do_short = 1;
    }

    if(do_short)
    {
        /* The long-range approximation isn't good enough */
        sum = 1.0;
        term = 1.0;

        i = 0;
        do {
            i++;
            term *= t2 / (2*m + 2*i + 1);
            test = sum;
            sum += term;
            test = sum - test;
        } while(test != 0.0);

        F[m] = sum * et / (2*m+1);
    }

    /* Now do downwards recursion */
    for(i = m - 1; i >= 0; i--)
        F[i] = (t2 * F[i + 1] + et) / (2 * i + 1);
}


static void mirp_boys_double_zero(double *F, int m)
{
    for(int i = 0; i <= m; i++)
        F[i] = 1.0 / (2.0 * (double)i + 1);
}


void mirp_boys_double(double *F, int m, double t)
{
    assert(m >= 0);
    assert(working_prec > 0);

    if(t == 0.0)
        mirp_boys_double_zero(F, m);
    else    
        mirp_boys_double_nonzero(F, m, t);
}

