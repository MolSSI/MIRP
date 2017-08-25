/*! \file
 *
 * \brief Some miscellaneous mathematical functions in double
 *        precision and interval arithmetic
 */

#include "mirp/math.h"
#include <assert.h>


slong mirp_min_accuracy_bits(arb_srcptr v, size_t n)
{
    if(n == 0)
        return 0;

    slong min = arb_rel_accuracy_bits(v+0);

    for(size_t i = 1; i < n; i++)
        min = MIN(min, arb_rel_accuracy_bits(v+i));

    return min;
}


double mirp_factorial_d(int n)
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


double mirp_factorial2_d(int n)
{
    assert(n >= -1);

    if(n >= -1 && n <= 1)
        return 1.0;
    else
    {
        double r = 1.0;
        for(int i = n; i > 0; i-=2)
            r *= i;
        return r;
    }
}


void mirp_factorial2(arb_t output, long int n, slong working_prec)
{
    /* arblib has a double factorial function, but only for
     * positive values */

    if(n >= -1 && n <= 1)
        arb_set_ui(output, 1);
    else
        arb_doublefac_ui(output, (unsigned long)n, working_prec);
}


double mirp_binomial_d(int n, int k)
{
    assert(n >= 0);
    assert(k >= 0);
    assert(k <= n);

    return mirp_factorial_d(n)/(mirp_factorial_d(k) * mirp_factorial_d(n-k));
}


void mirp_binomial(arb_t output, long int n, long int k, slong working_prec)
{
    assert(n >= 0);
    assert(k >= 0);
    assert(k <= n);

    /* A temporary */
    arb_t f;
    arb_init(f);

    /* put k! in the output */
    arb_fac_ui(output, k, working_prec);

    /* next, (n-k)! in temporary, and multiply into output */
    arb_fac_ui(f, n-k, working_prec);
    arb_mul(output, output, f, working_prec);

    /* now n! in temporary, and divide */
    arb_fac_ui(f, n, working_prec);
    arb_div(output, f, output, working_prec);

    arb_clear(f);
}

