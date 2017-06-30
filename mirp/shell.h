/*! \file
 *
 * \brief Functions for manipulating shell structures
 */

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

void mirp_normalize_shell_double(int am, int nprim, int ngeneral,
                                 const double * alpha,
                                 const double * coeff,
                                 double * coeff_out);

void mirp_normalize_shell_interval(int am, int nprim, int ngeneral,
                                   const arb_t * alpha,
                                   const arb_t * coeff,
                                   arb_t * coeff_out,
                                   slong working_prec);

int mirp_iterate_gaussian(int lmn[3]);

#ifdef __cplusplus
}
#endif

