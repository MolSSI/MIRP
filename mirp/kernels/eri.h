/*! \file
 *
 * \brief Calculation of electron repulsion integrals
 */

#pragma once

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Computes a single cartesian electron repulsion integral
 *         (double precision)
 *
 * The \p lmn paramters are the exponents on x, y, and z, with the total
 * angular momentum being the sum of the three components. For example,
 * { 1, 2, 0 } is the \f$f_{xy^2}\f$ cartesian integral.
 *
 * \param [out] output 
 *              Resulting integral output
 * \param [in]  lmn1,lmn2,lmn3,lmn4
 *              Exponents of x, y, and z that signify angular momentum. Required
 *              to be 3 elements.
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four centers
 */
void mirp_eri_single_double(double * output,
                            const int * lmn1, const double * A, double alpha1,
                            const int * lmn2, const double * B, double alpha2,
                            const int * lmn3, const double * C, double alpha3,
                            const int * lmn4, const double * D, double alpha4);


/*! \brief Compute all cartesian components of a single primitive electron
 *         repulsion integral (double precision)
 *
 * The \p output buffer is expected to be able to hold all primitive integrals
 * (ie, it can hold ncart(am1) * ncart(am2) * ncart(am3) * ncart(am4) elements).
 *
 * \param [out] output 
 *              Resulting integral output
 * \param [in]  am1,am2,am3,am4
 *              Angular momentum of the four centers
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four centers
 */
void mirp_prim_eri_double(double * output,
                          int am1, const double * A, double alpha1,
                          int am2, const double * B, double alpha2,
                          int am3, const double * C, double alpha3,
                          int am4, const double * D, double alpha4);


/*! \brief Compute all cartesian components of the electron repulsion integral
           of contracted shell quartet (four center, double precision)
 *
 * The coefficients (\p coeff1, \p coeff2, \p coeff3, \p coeff4) are expected
 * to be unnormalized.
 *
 * The \p output buffer is expected to be able to hold all the resulting
 * integrals
 *
 * For the coefficient arrays, the index of the primitive/exponent is the fastest index.
 * (The index is `[g*nprim+p]` where `g` is the index of the general contraction and `p`
 * the index of the primitive).
 *
 *
 * \param [out] output
 *              Completed integral output
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
 *              (lengths \p nprim1, \p nprim2, \p nprim3, \p nprim4)
 * \param [in]  coeff1,coeff2,coeff3,coeff4
 *              Coefficients for each primitive gaussian on the four centers
 *              (\p nprim1 * \p ngeneral1, \p nprim2 * \p ngeneral2,
 *              \p nprim3 * \p ngeneral3, \p nprim4 * \p ngeneral4)
 */
void mirp_eri_double(double * output,
                     int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1,
                     int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2,
                     int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3,
                     int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4);


/*! \brief Computes a single cartesian electron repulsion integral
 *         (interval arithmetic)
 *
 * \copydetails mirp_eri_single_double
 * \param [in] working_prec The working precision (binary digits/bits) to use in the calculation
 */
void mirp_eri_single_interval(arb_t output,
                              const int * lmn1, const arb_t * A, const arb_t alpha1,
                              const int * lmn2, const arb_t * B, const arb_t alpha2,
                              const int * lmn3, const arb_t * C, const arb_t alpha3,
                              const int * lmn4, const arb_t * D, const arb_t alpha4,
                              slong working_prec);


/*! \brief Compute all cartesian components of a single primitive electron
 *         repulsion integral (interval arithmetic)
 *
 * \copydetails mirp_prim_eri_double
 * \param [in] working_prec The working precision (binary digits/bits) to use in the calculation
 */
void mirp_prim_eri_interval(arb_t * output,
                            int am1, const arb_t * A, const arb_t alpha1,
                            int am2, const arb_t * B, const arb_t alpha2,
                            int am3, const arb_t * C, const arb_t alpha3,
                            int am4, const arb_t * D, const arb_t alpha4,
                            slong working_prec);


/*! \brief Compute all cartesian components of the electron repulsion integral
           of contracted shell quartet (interval arithmetic)
 *
 * \copydetails mirp_eri_double
 * \param [in] working_prec The working precision (binary digits/bits) to use in the calculation
 */
void mirp_eri_interval(arb_t * output,
                       int am1, const arb_t * A, int nprim1, int ngeneral1, const arb_t * alpha1, const arb_t * coeff1,
                       int am2, const arb_t * B, int nprim2, int ngeneral2, const arb_t * alpha2, const arb_t * coeff2,
                       int am3, const arb_t * C, int nprim3, int ngeneral3, const arb_t * alpha3, const arb_t * coeff3,
                       int am4, const arb_t * D, int nprim4, int ngeneral4, const arb_t * alpha4, const arb_t * coeff4,
                       slong working_prec);

#ifdef __cplusplus
}
#endif

