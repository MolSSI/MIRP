/*! \file
 *
 * \brief Some miscellaneous mathematical functions in double
 *        precision and interval arithmetic
 */

#pragma once

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Calculates the minimum value of two variables */
#define MIN(a,b) (((a)<(b))?(a):(b))


/*! \brief Calculates the maximum value of two variables */
#define MAX(a,b) (((a)>(b))?(a):(b))


/*! \brief Calculates (-1) to a given power
 *
 * This replaces pow(-1, n)
 */
#define NEG1_POW(n) (((n)%2)?(-1):(1))


/*! \brief The value of the constant pi in double precision */
#define MIRP_PI 3.14159265358979324


/*! \brief The value of the sqrt(pi) in double precision */
#define MIRP_SQRT_PI 1.7724538509055160273 


/*! \brief The value of the pi**(3/2) in double precision */
#define MIRP_PI_32 5.5683279968317078453
    

/*! \brief log_10(2)
 * 
 * For conversion between binary precision and decimal precision (number of digits)
 * ndigits = log10(2) * precision (in bits)
 */
#define MIRP_LOG_10_2 0.3010299956639812


/*! \brief Tests if a value with no significant bits is zero to a given precision
 *
 * If a number does not have significant bits, that means that it is 0 +/- [error].
 * This function first checks if there are significant bits. If there isn't,
 * it then checks whether or not the error associated with
 * this zero value is small enough that it would be completely zero
 * when represented with the given (binary) precision.
 *
 * \note This function uses an internal safety factor, so the results may be
 *       unexpected if you are testing.
 *
 * \param [in] n The value to test
 * \param [in] prec The binary precision to test to
 * \return Nonzero if the test value has no signifcant bits and is zero
 *         to the given precision
 */
int mirp_test_zero_prec(const arb_t n, slong prec);


/*! \brief Calculates a factorial using double precision */
double mirp_factorial(int n);


/*! \brief Calculates a double factorial using double precision */
double mirp_double_factorial(int n);


/*! \brief Calculates a double factorial using interval arithmetic */
void mirp_double_factorial_interval(arb_t output, long int n, slong working_prec);


/*! \brief Calculates a binomial coefficient using double precision */
double mirp_binomial_coefficient(int n, int k);


/*! \brief Calculates a binomial coefficient using interval arithmetic */
void mirp_binomial_coefficient_interval(arb_t output, long int n, long int k, slong working_prec);

#ifdef __cplusplus
}
#endif

