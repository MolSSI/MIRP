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

/*! \brief Pointer to a function that computes a single cartesian integral
 *         (four-center, double precision)
 */
typedef void (*cb_integral4_single_d)(double *,
                                      const int *, const double *, double,
                                      const int *, const double *, double,
                                      const int *, const double *, double,
                                      const int *, const double *, double);


/*! \brief Pointer to a function that computes all cartesian integrals
 *         for a contracted shell quartet (four-center, double precision)
 */
typedef void (*cb_integral4_d)(double *,
                               int, const double *, int, int, const double *, const double *,
                               int, const double *, int, int, const double *, const double *,
                               int, const double *, int, int, const double *, const double *,
                               int, const double *, int, int, const double *, const double *);



/*******************************************
 * Interval arithmetic
 *******************************************/

/*! \brief Pointer to a function that computes a single cartesian integral
 *         (four-center, interval arithmetic)
 */
typedef void (*cb_integral4_single)(arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    slong);


/*! \brief Pointer to a function that computes a single cartesian integral
 *         to a specified accuracy (four-center)
 */
typedef int (*cb_integral4_single_target)(arb_t,
                                          const int *, arb_srcptr, const arb_t,
                                          const int *, arb_srcptr, const arb_t,
                                          const int *, arb_srcptr, const arb_t,
                                          const int *, arb_srcptr, const arb_t,
                                          slong);


/*! \brief Pointer to a function that computes a single cartesian integral
 *         to a specified accuracy (four-center, with string inputs)
 */
typedef void (*cb_integral4_single_target_str)(arb_t,
                                               const int *, const char **, const char *,
                                               const int *, const char **, const char *,
                                               const int *, const char **, const char *,
                                               const int *, const char **, const char *,
                                               slong);


/*! \brief Pointer to a function that computes a single cartesian integral
 *         to exact double precision (four-center)
 */
typedef void (*cb_integral4_single_exact)(double * integral,
                                          const int * lmn1, const double * A, double alpha1,
                                          const int * lmn2, const double * B, double alpha2,
                                          const int * lmn3, const double * C, double alpha3,
                                          const int * lmn4, const double * D, double alpha4);


/*! \brief Pointer to a function that computes all cartesian integrals
 *         for a contracted shell quartet (four-center, interval arithmetic)
 */
typedef void (*cb_integral4)(arb_ptr,
                             int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                             int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                             int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                             int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                             slong);


/*! \brief Pointer to a function that computes all cartesian integrals
 *         for a contracted shell quartet to a specified accuracy (four-center)
 */
typedef int (*cb_integral4_target)(arb_ptr,
                                   int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                                   int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                                   int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                                   int, arb_srcptr, int, int, arb_srcptr, arb_srcptr,
                                   slong);


/*! \brief Pointer to a function that computes all cartesian integrals
 *         for a contracted shell quartet to a specified accuracy (four-center, with string inputs)
 */
typedef void (*cb_integral4_target_str)(arb_ptr,
                                        int, const char **, int, int, const char **, const char **,
                                        int, const char **, int, int, const char **, const char **,
                                        int, const char **, int, int, const char **, const char **,
                                        int, const char **, int, int, const char **, const char **,
                                        slong);


/*! \brief Pointer to a function that computes all cartesian integrals
 *         for a contracted shell quartet to exact double precision (four-center)
 */
typedef int (*cb_integral4_exact)(double *,
                                  int, const double *, int, int, const double *, const double *,
                                  int, const double *, int, int, const double *, const double *,
                                  int, const double *, int, int, const double *, const double *,
                                  int, const double *, int, int, const double *, const double *);

#ifdef __cplusplus
}
#endif

