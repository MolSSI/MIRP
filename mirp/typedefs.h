/*! \file
 *
 * \brief Typedefs of common function pointers
 */

#pragma once

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************
 * Double precision
 *******************************************/

/*! \brief A callback for a function that computes a single cartesian integral
 *         (four center, double precision)
 */
typedef void (*cb_integral4_single_d)(double *,
                                      const int *, const double *, double,
                                      const int *, const double *, double,
                                      const int *, const double *, double,
                                      const int *, const double *, double);


/*! \brief A callback for a function that computes all cartesian integrals for an
 *         entire shell (four center, double precision)
 */
typedef void (*cb_integral4_d)(double *,
                               int, const double *, int, int, const double *, const double *,
                               int, const double *, int, int, const double *, const double *,
                               int, const double *, int, int, const double *, const double *,
                               int, const double *, int, int, const double *, const double *);



/*******************************************
 * Interval arithmetic
 *******************************************/

/*! \brief A callback for a function that computes a single cartesian integral
 *         (four center, interval arithmetic)
 */
typedef void (*cb_integral4_single)(arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    slong);


typedef int (*cb_integral4_single_target)(arb_t,
                                          const int *, arb_srcptr, const arb_t,
                                          const int *, arb_srcptr, const arb_t,
                                          const int *, arb_srcptr, const arb_t,
                                          const int *, arb_srcptr, const arb_t,
                                          slong);


typedef void (*cb_integral4_single_target_str)(arb_t,
                                               const int *, const char **, const char *,
                                               const int *, const char **, const char *,
                                               const int *, const char **, const char *,
                                               const int *, const char **, const char *,
                                               slong);

typedef void (*cb_integral4_single_exact)(double * integral,
                                          const int * lmn1, const double * A, double alpha1,
                                          const int * lmn2, const double * B, double alpha2,
                                          const int * lmn3, const double * C, double alpha3,
                                          const int * lmn4, const double * D, double alpha4);


/*! \brief A callback for a function that computes all cartesian integrals for an
 *         entire shell (four center, interval arithmetic)
 */
typedef void (*cb_integral4)(arb_ptr,
                             int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                             int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                             int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                             int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                             slong);

#ifdef __cplusplus
}
#endif

