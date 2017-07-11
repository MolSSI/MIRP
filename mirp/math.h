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


/*! \brief The value of the pi**(3/2) in double precision */
#define MIRP_PI_32 5.5683279968317078453
    

/*! \brief log_10(2)
 * 
 * For conversion between binary precision and decimal precision (number of digits)
 * ndigits = log10(2) * precision (in bits)
 */
#define MIRP_LOG_10_2 0.3010299956639812


/*! \brief Number of cartesian functions for a given angular momentum */
#define MIRP_NCART(am) ((((am)+1)*((am)+2))/2)


/*! \brief Calculates a factorial using double precision */
double mirp_factorial(int n);


/*! \brief Calculates a double factorial using double precision */
double mirp_double_factorial(int n);


/*! \brief Calculates a double factorial using arbitrary precision */
void mirp_double_factorial_mp(mpfr_t result, long int n);


/*! \brief Calculates a double factorial using interval arithmetic */
void mirp_double_factorial_interval(arb_t result, long int n, slong working_prec);


/*! \brief Calculates a binomial coefficient using double precision */
double mirp_binomial_coefficient(int n, int k);


/*! \brief Calculates a binomial coefficient using arbitrary precision */
void mirp_binomial_coefficient_mp(mpfr_t result, long int n, long int k);


/*! \brief Calculates a binomial coefficient using interval arithmetic */
void mirp_binomial_coefficient_interval(arb_t result, long int n, long int k, slong working_prec);

#ifdef __cplusplus
}
#endif

