/*! \file
 *
 * \brief Functions related to gaussians and shells
 */

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Number of cartesian functions for a given angular momentum */
#define MIRP_NCART(am) ((((am)+1)*((am)+2))/2)

/*! \brief Number of cartesian functions for 2 shells */
#define MIRP_NCART2(am1, am2) \
        (MIRP_NCART((am1)) * MIRP_NCART((am2)))

/*! \brief Number of cartesian functions for 4 shells */
#define MIRP_NCART4(am1, am2, am3, am4) \
        (MIRP_NCART2((am1),(am2)) * MIRP_NCART2((am3),(am4)))

/*! \brief Number of cartesian functions for an lmn triplet */
#define MIRP_NCART_LMN(lmn) \
        (MIRP_NCART((lmn[0])+(lmn[1])+(lmn[2])))

/*! \brief Number of cartesian functions for 2 lmn triplets */
#define MIRP_NCART_LMN2(lmn1, lmn2) \
        (MIRP_NCART_LMN((lmn1)) * MIRP_NCART_LMN((lmn2)))

/*! \brief Number of cartesian functions for 4 lmn triplets */
#define MIRP_NCART_LMN4(lmn1, lmn2, lmn3, lmn4)  \
        (MIRP_NCART_LMN2((lmn1),(lmn2)) * MIRP_NCART_LMN2((lmn3),(lmn4)))


/*! \brief Find the next gaussian in the ordering
 *
 * Obtain the next l, m, and n parameters of a gaussian in the internal MIRP
 * ordering
 *
 * For example, if \p lmn = {2, 1, 0} is input, the result will be {2, 0, 1}.
 *
 * If the return value of this function is 0, the contents of \p lmn are not
 * defined.
 *
 * \param [inout] lmn The l, m, and n parameters of a gaussian basis function
 * \return Nonzero if the new \p lmn is a valid gaussian, 0 if it is not
 *         (i.e., we have iterated past the end of the complete set of
 *         gaussians)
 */
int mirp_iterate_gaussian(int * lmn);


/*! \brief Create all lmn combinations for a given angular momentum
 *
 * This fills in the \p lmn parameter will all combinations of l,m, and n
 * that are valid for the given angular momentum. They will be in order.
 *
 * \warning The \p lmn parameter must be allocated and large enough
 *          to hold 3*number of cartesian components for the given \p am
 *
 * \param [in] am   The angular momentum
 * \param [out] lmn Place to put all the lmn combinations
 */
void mirp_gaussian_fill_lmn(int am, int * lmn);


/*! \brief Normalize a shell (double precision)
 *
 * This function normalizes the contraction coefficients of the shell.
 *
 * \param [in] am         The angular momentum of the shell (0 = s, 1 = p, etc)
 * \param [in] nprim      Number of primitives in the shell
 * \param [in] ngeneral   Number of general contractions in the shell
 * \param [in] alpha      The exponents of the shell (length \p nprim)
 * \param [in] coeff      The (unnormalized) contraction coefficients
 *                        (length \p nprim * \p ngeneral)
 * \param [out] coeff_out Normalized contraction coefficients
 *                        (length \p nprim * \p ngeneral)
 */
void mirp_normalize_shell_d(int am, int nprim, int ngeneral,
                            const double * alpha,
                            const double * coeff,
                            double * coeff_out);


/*! \brief Normalize a shell (interval arithmetic)
 *
 * \copydetails mirp_normalize_shell_d
 * \param [in] working_prec The working precision (binary digits/bits) to use
 *                          in the calculation
 */
void mirp_normalize_shell(int am, int nprim, int ngeneral,
                          arb_srcptr alpha,
                          arb_srcptr coeff,
                          arb_ptr coeff_out,
                          slong working_prec);


#ifdef __cplusplus
}
#endif

