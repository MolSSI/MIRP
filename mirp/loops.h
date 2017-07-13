/*! \file
 *
 * \brief Helpers for looping over cartesian and primitive integrals
 */

#pragma once

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************
 * Double precision
 *******************************************/

/*! \brief A callback for a function that computes a single cartesian integral (double precision) */
typedef void (*cb_single4_double)(double *,
                                  const int *, const double *, double,
                                  const int *, const double *, double,
                                  const int *, const double *, double,
                                  const int *, const double *, double);


/*! \brief A callback for a function that computes all cartesian integrals for a primitive quartet (double precision) */
typedef void (*cb_prim4_double)(double *,
                                int, const double *, double,
                                int, const double *, double,
                                int, const double *, double,
                                int, const double *, double);


/*! \brief Compute all cartesian components of a single primitive integral (double precision)
 *
 * The \p result buffer is expected to be able to hold all primitive integrals (ie, it
 * can hold ncart(am1) * ncart(am2) * ncart(am3) * ncart(am4) elements).
 *
 * \param [out] result 
 *              Resulting integral output
 * \param [in]  am1,am2,am3,am4
 *              Angular momentum of the four centers
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four centers
 * \param [in]  cb
 *              Callback that calculates a single cartesian component of a primitive integral
 */
void mirp_cartloop4_double(double * result,
                           int am1, const double * A, double alpha1,
                           int am2, const double * B, double alpha2,
                           int am3, const double * C, double alpha3,
                           int am4, const double * D, double alpha4,
                           cb_single4_double cb);


/*! \brief Compute all cartesian components of a contracted shell (double precision)
 *
 * The coefficients (\p coeff1, \p coeff2, \p coeff3, \p coeff4) are expected to be unnormalized.
 *
 * The \p result buffer is expected to be able to hold all the resulting integrals
 *
 * \param [out] result
 *              Resulting integral output
 * \param [in]  am1,am2,am3,am4
 *              Angular momentum of the four centers
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  nprim1,nprim2,nprim3,nprim4
 *              Number of primitives on each center
 * \param [in]  ngeneral1,ngeneral2,ngeneral3,ngeneral4
 *              Number of general contractions on each center
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the primitive gaussians on the four centers
 *              (lengths \p nprim1, \p nprim2, \p nprim3, \p nprim4, respectively)
 * \param [in]  coeff1,coeff2,coeff3,coeff4
 *              Coefficients for each primitive gaussian on the four centers (lengths
 *              \p nprim1 * \p ngeneral1, \p nprim2 * \p ngeneral2,
 *              \p nprim3 * \p ngeneral3, \p nprim4 * \p ngeneral4, respectively)
 * \param [in]  cb
 *              Callback that calculates all cartesian components of a single primitive integral
 */
void mirp_loop4_double(double * result,
                       int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1,
                       int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2,
                       int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3,
                       int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4,
                       cb_prim4_double cb);



/*******************************************
 * Interval arithmetic
 *******************************************/

/*! \brief A callback for a function that computes a single cartesian integral (arbitrary precision) */
typedef void (*cb_single4_interval)(arb_t,
                                    const int *, const arb_t *, const arb_t,
                                    const int *, const arb_t *, const arb_t,
                                    const int *, const arb_t *, const arb_t,
                                    const int *, const arb_t *, const arb_t,
                                    slong);


/*! \brief A callback for a function that computes all cartesian integrals for a primitive quartet (arbitrary precision) */
typedef void (*cb_prim4_interval)(arb_t *,
                                  int, const arb_t *, const arb_t,
                                  int, const arb_t *, const arb_t,
                                  int, const arb_t *, const arb_t,
                                  int, const arb_t *, const arb_t,
                                  slong);


/*! \brief Compute all cartesian components of a single primitive integral (interval arithmetic)
 *
 * \copydetails mirp_cartloop4_double
 * \param [in] working_prec The working precision to use in the calculation
 */
void mirp_cartloop4_interval(arb_t * result,
                             int am1, const arb_t * A, const arb_t alpha1,
                             int am2, const arb_t * B, const arb_t alpha2,
                             int am3, const arb_t * C, const arb_t alpha3,
                             int am4, const arb_t * D, const arb_t alpha4,
                             slong working_prec, cb_single4_interval cb);


/*! \brief Compute all cartesian components of a contracted shell (interval arithmetic)
 *
 * \copydetails mirp_loop4_double
 * \param [in] working_prec The working precision to use in the calculation
 */
void mirp_loop4_interval(arb_t * result,
                         int am1, const arb_t * A, int nprim1, int ngeneral1, const arb_t * alpha1, const arb_t * coeff1,
                         int am2, const arb_t * B, int nprim2, int ngeneral2, const arb_t * alpha2, const arb_t * coeff2,
                         int am3, const arb_t * C, int nprim3, int ngeneral3, const arb_t * alpha3, const arb_t * coeff3,
                         int am4, const arb_t * D, int nprim4, int ngeneral4, const arb_t * alpha4, const arb_t * coeff4,
                         slong working_prec, cb_prim4_interval cb);


#ifdef __cplusplus
}
#endif

