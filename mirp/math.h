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

/*! \brief Smallest subnormal double-precision value */
#define MIRP_DBL_TRUE_MIN 4.9406564584124654e-324


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
 * For conversion between binary precision and decimal precision
 * (number of decimal digits)
 *
 * ndigits = log10(2) * precision (in bits)
 */
#define MIRP_LOG_10_2 0.3010299956639812


/*! \brief Finds the least number of accuracy bits in a vector
 *
 * \param [in] v The vector to search the minimum accuracy of
 * \param [in] n The length of the vector
 * \return The minimum accuracy found of all elements of the vector
 */
slong mirp_min_accuracy_bits(arb_srcptr v, size_t n);


/*! \brief Calculates b^e with e being a signed integer */
void mirp_pow_si(arb_t output, const arb_t b, long e, slong prec);


/*! \brief Calculates a factorial using interval arithmetic */
void mirp_factorial(arb_t output, long n);


/*! \brief Calculates a factorial using double precision */
double mirp_factorial_d(int n);


/*! \brief Calculates a double factorial using double precision */
double mirp_factorial2_d(int n);


/*! \brief Calculates a double factorial using interval arithmetic */
void mirp_factorial2(arb_t output, long n);


/*! \brief Calculates a binomial coefficient using double precision */
double mirp_binomial_d(int n, int k);


/*! \brief Calculates a binomial coefficient using interval arithmetic */
void mirp_binomial(arb_t output, long n, long k);


#ifdef __cplusplus
}
#endif

