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


void mirp_binomial_coefficient_interval(arb_t result, long int n, long int k, slong working_prec)
{
    assert(n >= 0);
    assert(k >= 0);
    assert(k <= n);

    /* A temporary */
    arb_t f;
    arb_init(f);

    /* put k! in the result */
    arb_fac_ui(result, k, working_prec);

    /* next, (n-k)! in temporary, and multiply into result */
    arb_fac_ui(f, n-k, working_prec);
    arb_mul(result, result, f, working_prec);

    /* now n! in temporary, and divide */
    arb_fac_ui(f, n, working_prec);
    arb_div(result, f, result, working_prec);

    arb_clear(f);
}

