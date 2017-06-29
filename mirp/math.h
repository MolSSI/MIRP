/*! \file
 *
 * \brief Some miscellaneous mathematical functions in double
 *        precision and arbitrary precision
 */

#pragma once

#include <gmp.h>
#include <mpfr.h>
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
    
/* \brief log_10(2)
 * 
 * For conversion between binary precision and decimal precision (number of digits)
 * ndigits = log10(2) * precision (in bits)
 */
#define MIRP_LOG_10_2 0.3010299956639812


/* \brief Number of cartesian functions for a given angular momentum */
#define MIRP_NCART(am) ((((am)+1)*((am)+2))/2)


/* \brief Calculates a factorial in double precision
 *
 * A factorial is defined as
 *
 * n! = n*(n-1)*(n-2)*...*3*2*1
 * 0! = 1
 * 1! = 1
 *
 * \param [in] n The value of which to calculate the factorial
 * \return The resulting factorial
 */
double mirp_factorial(int n);


/*! \brief Calculates a double factorial in double precision
 *
 * A double factorial is defined as
 *
 * n!! = n*(n-2)*(n-4)*...*6*4*2  for n > 0 and even n
 * n!! = n*(n-2)*(n-4)*...*5*3*1  for n > 0 and odd n
 *    1!! = 1
 *    0!! = 1
 * (-1)!! = 1
 *
 * \warning Double factorial is defined for negative numbers. However, this
 *          function will only calculate n!! for n >= -1
 *
 * \param [in] n The value of which to calculate the double factorial
 * \return The resulting double factorial
 */
double mirp_double_factorial(int n);


/*! \brief Calculates a double factorial in arbitrary precision
 *
 * A double factorial is defined as
 *
 * n!! = n*(n-2)*(n-4)*...*6*4*2  for n > 0 and even n
 * n!! = n*(n-2)*(n-4)*...*5*3*1  for n > 0 and odd n
 *    1!! = 1
 *    0!! = 1
 * (-1)!! = 1
 *
 * The working precision is taken to be the same as \p result
 *
 * \warning Double factorial is defined for negative numbers. However, this
 *          function will only calculate n!! for n >= -1
 *
 * \param [in] n The value of which to calculate the double factorial
 * \param [out] result The resulting double factorial
 */
void mirp_double_factorial_mp(mpfr_t result, long int n);


/*! \brief Calculates a double factorial with interval arithmetic 
 *
 * A double factorial is defined as
 *
 * n!! = n*(n-2)*(n-4)*...*6*4*2  for n > 0 and even n
 * n!! = n*(n-2)*(n-4)*...*5*3*1  for n > 0 and odd n
 *    1!! = 1
 *    0!! = 1
 * (-1)!! = 1
 *
 * The working precision is taken to be the same as \p result
 *
 * \warning Double factorial is defined for negative numbers. However, this
 *          function will only calculate n!! for n >= -1
 *
 * \param [in] n The value of which to calculate the double factorial
 * \param [out] result The resulting double factorial
 */
void mirp_double_factorial_interval(arb_t result, long int n, slong working_prec);


/*! \brief Calculates a binomial coefficient in double precision
 *
 * A binomial coefficient is defined as
 *
 * (n,k) = n! / (k! * (n-k)!)
 *
 * for n,k >= 0 and k <= n
 *
 * \param [in] n The coefficient n
 * \param [in] k The coefficient k
 * \return The resulting binomial coefficient (n,k)
 */
double mirp_binomial_coefficient(int n, int k);


/*! \brief Calculates a binomial coefficient in arbitrary precision
 *
 * A binomial coefficient is defined as
 *
 * (n,k) = n! / (k! * (n-k)!)
 *
 * for n,k >= 0 and k <= n
 *
 * The working precision is taken to be the same as \p result
 *
 * \param [in] n The coefficient n
 * \param [in] k The coefficient k
 * \param [out] result The resulting binomial coefficient (n,k)
 */
void mirp_binomial_coefficient_mp(mpfr_t result, long int n, long int k);


/*! \brief Calculates a binomial coefficient with interval arithmetic
 *
 * A binomial coefficient is defined as
 *
 * (n,k) = n! / (k! * (n-k)!)
 *
 * for n,k >= 0 and k <= n
 *
 * The working precision is taken to be the same as \p result
 *
 * \param [in] n The coefficient n
 * \param [in] k The coefficient k
 * \param [out] result The resulting binomial coefficient (n,k)
 */
void mirp_binomial_coefficient_interval(arb_t result, long int n, long int k, slong working_prec);

#ifdef __cplusplus
}
#endif

