/*! \file
 *
 * \brief Some miscellaneous mathematical functions in
 *        interval arithmetic
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

void mirp_pow_si(arb_t output, const arb_t b, long e, slong prec)
{
    if(e >= 0)
        arb_pow_ui(output, b, (unsigned long)e, prec);
    else
    {
        arb_pow_ui(output, b, (unsigned long)(-e), prec);
        arb_ui_div(output, 1, output, prec);
    }
}


void mirp_factorial(arb_t output, long n)
{
    assert(n >= 0);

    if(n == 0 || n == 1)
        arb_one(output);
    else
    {
        fmpz_t tmp;
        fmpz_init_set_ui(tmp, 1ul);
        fmpz_fac_ui(tmp, (unsigned long)n); 
        arb_set_fmpz(output, tmp);
        fmpz_clear(tmp);
    }
}


void mirp_factorial2(arb_t output, long n)
{
    /* arblib has a double factorial function, but only for
     * positive values */

    assert(n >= -1);

    if(n >= -1 && n <= 1)
        arb_one(output);
    else
    {
        fmpz_t tmp;
        fmpz_init_set_ui(tmp, 1ul);

        for(long i = n; i > 1; i -= 2)
            fmpz_mul_ui(tmp, tmp, (unsigned long)i);

        arb_set_fmpz(output, tmp);
        fmpz_clear(tmp);
    }

}


void mirp_binomial(arb_t output, long n, long k)
{
    assert(n >= 0);
    assert(k >= 0);
    assert(k <= n);

    fmpz_t tmp;
    fmpz_init(tmp);
    fmpz_bin_uiui(tmp, (unsigned long)n, (unsigned long)k);
    arb_set_fmpz(output, tmp);
    fmpz_clear(tmp);
}

