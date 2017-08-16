/*! \file
 *
 * \brief Some miscellaneous mathematical functions in double
 *        precision and arbitrary precision
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


#if 0
int mirp_test_zero_prec(const arb_t n, slong prec)
{
    /* This function is meant to be used if n does not
       have any significant bits (ie, is 0 +/- some number)
     */
    if(arb_rel_accuracy_bits(n) > 0)
        return 0;

    /* This is needed to test if the result is zero to the number of digits
       we want (plus a safety factor, of course)
    
       The safety factor is 16 decimal digits / 53 binary digits
    */
    arb_t zero, zero_err;
    arb_init(zero);
    arb_init(zero_err);
    arb_zero(zero);
    arb_ui_pow_ui(zero_err, 10, 16, prec+53);
    arb_mul_2exp_si(zero_err, zero_err, prec+53);
    arb_ui_div(zero_err, 1, zero_err, prec+53);
    arb_add_error(zero, zero_err);

    int is_zero = arb_contains(zero, n);

    arb_clear(zero);
    arb_clear(zero_err);

    return is_zero;
}
#endif

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


void mirp_factorial2(arb_t output, long int n, slong working_prec)
{
    /* arblib has a double factorial function, but only for
     * positive values
     */
    
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

