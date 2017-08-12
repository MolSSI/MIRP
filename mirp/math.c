/*! \file
 *
 * \brief Some miscellaneous mathematical functions in double
 *        precision and arbitrary precision
 */

#include "mirp/math.h"
#include <assert.h>


/*! \brief Finds the minimum relative accuracy bit of a vector
 *
 * \param [in] v Vector to check
 * \param [in] n Length of the vector
 * \return The minimum value of arb_rel_accuracy_bits of the entire vector
 */
slong mirp_min_rel_accuracy_bits(arb_srcptr v, size_t n)
{
    if(n == 0)
        return -1;

    slong min = arb_rel_accuracy_bits(v+0);

    for(size_t i = 1; i < n; i++)
        min = MIN(min, arb_rel_accuracy_bits(v+i));

    return min;
}


/*! \brief Checks all elements of a vector for sufficient accuracy
 *
 * \param [in] v Vector to check
 * \param [in] n Length of the vector
 * \param [in] target_prec Precision required of all the elements
 * \return Nonzero if all elements are of sufficient accuracy, zero otherwise
 */
int mirp_all_sufficient_accuracy(arb_srcptr v, size_t n, slong target_prec)
{
    slong min_bits = mirp_min_rel_accuracy_bits(v, n);
    return min_bits >= target_prec;
}


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


void mirp_double_factorial_interval(arb_t output, long int n, slong working_prec)
{
    /* arblib has a double factorial function, but only for
     * positive values
     */
    
    if(n >= -1 && n <= 1)
        arb_set_ui(output, 1);
    else
        arb_doublefac_ui(output, (unsigned long)n, working_prec);
}


double mirp_binomial_coefficient(int n, int k)
{
    assert(n >= 0);
    assert(k >= 0);
    assert(k <= n);

    return mirp_factorial(n)/(mirp_factorial(k) * mirp_factorial(n-k));
}


void mirp_binomial_coefficient_interval(arb_t output, long int n, long int k, slong working_prec)
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

