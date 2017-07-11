/*! \file
 *
 * \brief Functions for manipulating shell structures
 */

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Normalize a shell (double precision)
 *
 * Obtain the next l, m, and n parameters of a gaussian in the internal MIRP ordering
 *
 * For example, if \p lmn = {2, 1, 0} is input, the result will be {2, 0, 1}.
 *
 * If the return value of this function is 0, the contents of \p lmn are not defined.
 *
 * \param [inout] lmn The l, m, and n parameters of a gaussian basis function
 * \return 1 if the new \p lmn is a valid gaussian, 0 if it is not
 *         (i.e., we have iterated past the end)
 */
int mirp_iterate_gaussian(int * lmn);


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
void mirp_normalize_shell_double(int am, int nprim, int ngeneral,
                                 const double * alpha,
                                 const double * coeff,
                                 double * coeff_out);


/*! \brief Normalize a shell (interval arithmetic)
 *
 * \copydetails mirp_normalize_shell_double
 * \param [in] working_prec The working precision to use in the calculation
 */
void mirp_normalize_shell_interval(int am, int nprim, int ngeneral,
                                   const arb_t * alpha,
                                   const arb_t * coeff,
                                   arb_t * coeff_out,
                                   slong working_prec);


#ifdef __cplusplus
}
#endif

