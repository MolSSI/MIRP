/*! \file
 *
 * \brief Some miscellaneous mathematical functions in double
 *        precision and arbitrary precision
 */

#include "mirp/math.h"
#include <assert.h>

double mirp_factorial(int n)
{
    assert(n >= 0);

    if(n == 0)
        return 1.0;
    if(n == 1)
        return 1.0;
    else
    {
        double r = 1.0;
        for(int i = n; i > 0; i--)
            r *= i;
        return r;
    }
}


double mirp_double_factorial(int n)
{
    assert(n >= -1);

    if(n == -1)
        return 1.0;
    if(n == 0)
        return 1.0;
    if(n == 1)
        return 1.0;
    else
    {
        double r = 1.0;
        for(int i = n; i > 0; i-=2)
            r *= i;
        return r;
    }
}


void mirp_double_factorial_mp(mpfr_t result, long int n)
{
    assert(n >= -1);

    /* Set the starting value */
    mpfr_set_ui(result, 1, MPFR_RNDN);

    /* Handle -1, 0, and 1 as special cases
     *
     * By definition, these all result in 1. Since result
     * was set to 1 above, we can just return
     */
    if(n >= -1 && n <= 1)
        return;

    for(long i = n; i > 0; i-=2)
        mpfr_mul_si(result, result, i, MPFR_RNDN);
}

void mirp_double_factorial_interval(arb_t result, long int n, slong working_prec)
{
    /* arblib has a double factorial function, but only for
     * positive values
     */
    
    if(n >= -1 && n <= 1)
        arb_set_ui(result, 1);
    else
        arb_doublefac_ui(result, (unsigned long)n, working_prec);
}


double mirp_binomial_coefficient(int n, int k)
{
    assert(n >= 0);
    assert(k >= 0);
    assert(k <= n);

    return mirp_factorial(n)/(mirp_factorial(k) * mirp_factorial(n-k));
}


void mirp_binomial_coefficient_mp(mpfr_t result, long int n, long int k)
{
    assert(n >= 0);
    assert(k >= 0);
    assert(k <= n);

    /* how many digits of precision are we using */
    mpfr_prec_t prec = mpfr_get_prec(result);

    /* A temporary */
    mpfr_t f;
    mpfr_init2(f, prec);

    /* put k! in the result */
    mpfr_fac_ui(result, k, MPFR_RNDN);

    /* next, (n-k)! in temporary, and multiply into result */
    mpfr_fac_ui(f, n-k, MPFR_RNDN);
    mpfr_mul(result, result, f, MPFR_RNDN);

    /* now n! in temporary, and divide */
    mpfr_fac_ui(f, n, MPFR_RNDN);
    mpfr_div(result, f, result, MPFR_RNDN);

    mpfr_clear(f);
}

