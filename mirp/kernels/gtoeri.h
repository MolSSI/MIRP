/*! \file
 *
 * \brief Kernel for electron repulsion integrals of gaussian orbitals
 */

#pragma once

#include <arb.h>
#include "mirp/kernels/integral4_wrappers.h"

#ifdef __cplusplus
extern "C" {
#endif



/*! \brief Computes a single cartesian GTO electron repulsion integral
 *         (interval arithmetic)
 *
 * \copydetails mirp_gtoeri_single_exact
 *
 * \param [in]  working_prec
 *              The working precision (binary digits/bits) to use
 *              in the calculation
 */
void mirp_gtoeri_single(arb_t integral,
                        const int * lmn1, arb_srcptr A, const arb_t alpha1,
                        const int * lmn2, arb_srcptr B, const arb_t alpha2,
                        const int * lmn3, arb_srcptr C, const arb_t alpha3,
                        const int * lmn4, arb_srcptr D, const arb_t alpha4,
                        slong working_prec);


/*******************
 * Wrappings
 *******************/

/*! \brief Compute a single GTO electron repulsion integral for a primitive quartet
 *         (interval arithmetic)
 *
 * \copydetails mirp_gtoeri_single
 */
MIRP_WRAP_SINGLE4_STR(gtoeri)


/*! \brief Compute a single GTO electron repulsion integral for a primitive quartet
 *         (exact double precision)
 *
 * The \p lmn parameters are the exponents on x, y, and z, with the total
 * angular momentum being the sum of the three components. For example,
 * { 1, 2, 0 } is the \f$f_{xy^2}\f$ cartesian integral.
 *
 * \param [out] integral
 *              Resulting integral integral
 * \param [in]  lmn1,lmn2,lmn3,lmn4
 *              Exponents of x, y, and z that signify angular momentum. Required
 *              to be 3 elements.
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four centers
 */
MIRP_WRAP_SINGLE4_EXACT(gtoeri)


/*! \brief Compute GTO electron repulsion integrals for a contracted
 *         shell quartet (interval arithmetic)
 *
 * \copydetails mirp_gtoeri_exact
 * \param [in]  working_prec
 *              The working precision (binary digits/bits) to use
 *              in the calculation
 */
MIRP_WRAP_SHELL4(gtoeri)


/*! \brief Compute GTO electron repulsion integrals for a contracted
 *         shell quartet (string inputs)
 *
 * \copydetails mirp_gtoeri
 */
MIRP_WRAP_SHELL4_STR(gtoeri)


/*! \brief Compute GTO electron repulsion integrals for a contracted
 *         shell quartet (exact double precision)
 *
 * \param [out] integrals
 *              Output for the computed integral
 * \param [in]  am1,am2,am3,am4
 *              Angular momentum for the four centers
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  nprim1,nprim2,nprim3,nprim4
 *              Number of primitive gaussians for each shell
 * \param [in]  ngen1,ngen2,ngen3,ngen4
 *              Number of general contractions for each shell
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the primitive gaussians on the four centers
 *              (of lengths \p nprim1, \p nprim2, \p nprim3, \p nprim4 respectively)
 * \param [in]  coeff1,coeff2,coeff3,coeff4
 *              Coefficients for all primitives and for all general contractions
 *              for each shell (of lengths \p nprim1 * \p ngen1, \p nprim2 * \p ngen2,
 *              \p nprim3 * \p ngen3, \p nprim4 * \p ngen4 respectively)
 */
MIRP_WRAP_SHELL4_EXACT(gtoeri)




#ifdef __cplusplus
}
#endif

